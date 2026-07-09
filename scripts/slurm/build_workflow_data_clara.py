import argparse
import json
import collections
import subprocess
from pathlib import Path
import os
from typing import List

# ---------------------------------------------------------------------------
# build_workflow_data_clara.py
#
# SWIF2 workflow generator for RECONSTRUCTION OF REAL (decoded) DATA using
# run-clara (multi-threaded, multi-file) instead of single-file recon-util.
#
# Each job processes a BATCH of FILES_PER_JOB decoded files together and
# reconstructs them with run-clara on THREADS threads. Per job:
#
#   1. SWIF2 stages the batch's decoded files from /mss to the worker node.
#   2. A single run-clara call reconstructs all staged files at once:
#         run-clara -y <yaml> -t THREADS -o <recon dir> -p rec_ <file>...
#      producing <recon dir>/rec_<original-basename> per input.
#
# run-clara needs CLARA_HOME pointed at the CLARA build paired with this
# coatjava; it is exported in the job command and passed via -c.
#
# Run this ON THE FARM (ifarm) where /mss is visible. For an off-farm preview,
# pass --data-dir pointing at a directory of placeholder files.
# ---------------------------------------------------------------------------

__JSONFORMAT = {"indent": 2, "separators": (",", ": ")}

RUN: str = "022083"

# Re-cook version suffix: bump this (e.g. "_2", "_3") to give every workflow a fresh
# name + /volatile output dir, so re-cooks don't collide with a previous swif2
# workflow of the same name or clobber its recon output. "" = no suffix.
RECOOK_SUFFIX: str = ""

# Input recon files to RE-COOK unconstrained. MUST still carry AHDC::wf (waveforms):
# fully-cooked DSTs without wf have AHDC recon silently skipped. The rich_060526 run
# dirs carry wf (validated). Change the run via RUN above (files are named rec_clas_<RUN>...).
DATA_DIR: str = f"/volatile/clas12/rg-l/production/rich_060526/calib/recon/{RUN}"

# Base /volatile area; the workflow output dir lives below this.
OUTPUT_BASE: str = "/volatile/clas12/touchte"

# Install base: the rgl-beamspot-rs tree holding the locally-built coatjava
# (with the ALERT beamConstraint flag) + its paired CLARA build and recon YAML.
INSTALL_BASE: str = "/work/clas12/users/touchte"

# Recon config: the new config.yaml (registers ALERTEngine, sets ALERT.beamConstraint="false").
YAML_PATH: str = f"{INSTALL_BASE}/data/simu/config.yaml"

# run-clara (multi-threaded) + the matching CLARA build for this coatjava.
COATJAVA_BIN: str = f"{INSTALL_BASE}/sync-coatjava/coatjava/bin"
RUN_CLARA: str = f"{COATJAVA_BIN}/run-clara"
CLARA_HOME: str = f"{INSTALL_BASE}/coatjava/clara"

# Batching / threading for run-clara.
FILES_PER_JOB: int = 2
THREADS: int = 24
RECON_PREFIX: str = "rec_"

DEBUG_COUNT: int = 2          # --debug: number of files for a quick validation
DEBUG_EVENTS: int = 1000      # --debug: events per file (run-clara -n); -1 = all
DEBUG_PARTITION: str = "priority"
DEBUG_TIME_SECS: str = "3600" # --debug: short walltime; the priority partition rejects long requests

# ---------------------------------------------------------------------------
# Single clean module-setup prefix (same as the other generators).
# ---------------------------------------------------------------------------
base_command: str = (
    "echo $SHELL; "
    "source /etc/profile.d/modules.sh; "
    "module use /scigroup/cvmfs/hallb/clas12/sw/modulefiles; "
    "module use /scigroup/cvmfs/geant4/modules; "
    "module purge; "
    "module load clas12; "
    'echo "Modules in use:" ; '
    "module list; "
)

# run-clara at 16 threads needs cores + RAM well above the single-file recon job.
recon_resources: dict[str, str] = {
    "constraint": "el9",
    "account": "clas12",
    "partition": "production",
    # NB: no "track": "osg" — inputs are site-local /volatile paths (not mss:), which
    # OSG cannot see, so off-site launch fails (SITE_LAUNCH_FAIL). Run on the JLab farm.
    "shell": "/bin/bash",
    "cpu_cores": str(THREADS),
    "ram_bytes": "34359738368",   # 32 GB
    "time_secs": "82800",         # 23 h — under the 24 h partition cap that rejected 86400
    "disk_bytes": "53687091200",  # 50 GB (5 staged inputs + recon output)
}


def remote_uri(path: str) -> str:
    """SWIF2 input URI for a staged file.

    Tape files (under /mss) need the `mss:` scheme; files already on disk
    (e.g. /volatile, /cache) are referenced by their bare path.
    """
    return f"mss:{path}" if path.startswith("/mss") else path


def with_partition(resources: dict[str, str], partition: str) -> dict[str, str]:
    """Return a copy of a resource dict with the partition overridden."""
    updated = dict(resources)
    updated["partition"] = partition
    return updated


def workflow_name_for() -> str:
    """Workflow / output-dir name for this run (with the re-cook suffix)."""
    return f"rgl_recon_unc_data_{RUN}{RECOOK_SUFFIX}"


def read_data_files(data_dir: str) -> List[str]:
    """List the decoded .hipo files for the run, sorted by chunk."""
    folder_path: Path = Path(data_dir)
    if folder_path.exists() and folder_path.is_dir():
        hipo_files: list[str] = [str(f) for f in folder_path.glob("*.hipo") if f.is_file()]
        hipo_files.sort()
        print(f"Found {len(hipo_files)} .hipo files in: {folder_path}")
        return hipo_files
    print(f"Folder does not exist or is not a directory: {folder_path}")
    return []


def chunk(items: List[str], size: int) -> List[List[str]]:
    """Split a list into consecutive batches of at most `size` items."""
    return [items[i:i + size] for i in range(0, len(items), size)]


def write_calibration_yaml(base_yaml: str, dest_yaml: str,
                           beam_x: float, beam_y: float) -> None:
    """Derive a recon YAML from `base_yaml` that adds an ALERT service block
    pinning the AHDC Kalman-filter beamline constraint to the assumed beam
    position (beam_x, beam_y) in mm, and write it to `dest_yaml`.

    The ALERT block goes under `configuration.services` (a sibling of AHDC), the
    block ALERTEngine reads via getEngineConfigString("beamX"/"beamY"). Inserted
    by a literal text splice (the repo does not depend on a YAML library); the
    base file is expected to be the unmodified config.yaml.
    """
    text = Path(base_yaml).read_text()
    if "\n    ALERT:" in text:
        raise SystemExit(
            f"{base_yaml} already contains an ALERT service block; point "
            f"--base-yaml at the clean config.yaml or edit beamX/beamY by hand.")
    anchor = '    AHDC:\n      Mode: "AI_GNN"\n'
    if anchor not in text:
        raise SystemExit(
            f"could not find the AHDC service block in {base_yaml}; insert the "
            f"ALERT beamX/beamY block by hand under configuration.services.")
    # Values MUST be quoted strings: CLAS12 ReconstructionEngine reads every engine
    # config value via JSONObject.getString (ALERTEngine then Double.parseDouble's it),
    # so an unquoted YAML number (beamX: 2.0) throws 'JSONObject["beamX"] not a string'.
    alert_block = (
        "    ALERT:\n"
        f'      beamX: "{beam_x}"\n'   # mm, assumed beam x (AHDC KF beamline constraint)
        f'      beamY: "{beam_y}"\n'   # mm, assumed beam y
    )
    Path(dest_yaml).write_text(text.replace(anchor, anchor + alert_block, 1))
    print(f"  wrote calibration YAML: {dest_yaml}  (ALERT beamX={beam_x}, beamY={beam_y} mm)")


def build_command(mss_files: List[str], output_recon: str, yaml_path: str,
                  n_events: int = -1) -> str:
    """run-clara the whole staged batch into `output_recon` with `yaml_path`.

    `output_recon` must be unique per job (run-clara writes a per-run
    filelist.txt there). `n_events` caps events per file via run-clara -n; -1
    (the default) omits the flag = all events (run-clara rejects an explicit
    "-n -1").
    """
    local_names = [Path(f).name for f in mss_files]

    nevents_opt = f"-n {n_events} " if n_events != -1 else ""
    inputs = " ".join(local_names)
    steps: list[str] = [
        f"mkdir -p {output_recon}",
        f"export CLARA_HOME={CLARA_HOME}",
        f"{RUN_CLARA} -y {yaml_path} -t {THREADS} {nevents_opt}-o {output_recon} "
        f"-p {RECON_PREFIX} -c {CLARA_HOME} {inputs}",
    ]
    return "; ".join(steps)


def create_recon_job(name: str, command: str, mss_files: List[str], total: int,
                     index: int, resources: dict[str, str]) -> collections.OrderedDict:
    """Build one batched run-clara job, staging all inputs to the worker node.

    Every file in the batch is declared as a SWIF2 input so it is staged to the
    node under its basename before the command runs. Tape inputs (/mss) use the
    `mss:` scheme; disk inputs (/volatile, /cache) use their bare path.
    """
    full_command = base_command + f'echo "Processing batch {index + 1}/{total}"; ' + command

    job = collections.OrderedDict(resources)
    job.update({
        "name": name,
        "phase": 0,
        "inputs": [{"local": Path(f).name, "remote": remote_uri(f)} for f in mss_files],
        "command": [full_command],
    })
    return job


def submit_workflow(workflow_name: str, json_path: Path) -> None:
    """Import and run a generated workflow JSON via swif2.

    Runs `swif2 import -file <json>` then `swif2 run -workflow <name>`. Must be
    executed on a node where the swif2 client is available (ifarm).
    """
    import_cmd = ["swif2", "import", "-file", str(json_path)]
    run_cmd = ["swif2", "run", "-workflow", workflow_name]
    for cmd in (import_cmd, run_cmd):
        print(f"  $ {' '.join(cmd)}")
        result = subprocess.run(cmd)
        if result.returncode != 0:
            print(f"  ERROR: command failed (exit {result.returncode}); "
                  f"skipping remaining submit steps for {workflow_name}")
            return


def main():
    parser = argparse.ArgumentParser(
        description="Generate SWIF2 workflow JSON: run-clara recon of real decoded data")
    parser.add_argument("--debug", action="store_true",
                        help=f"Small validation workflow ({DEBUG_COUNT} files) on the {DEBUG_PARTITION} partition")
    parser.add_argument("--debug-count", type=int, default=DEBUG_COUNT,
                        help=f"Number of files in --debug mode (default: {DEBUG_COUNT})")
    parser.add_argument("--debug-events", type=int, default=DEBUG_EVENTS,
                        help=f"Events per file in --debug mode (default: {DEBUG_EVENTS})")
    parser.add_argument("--files-per-job", type=int, default=FILES_PER_JOB,
                        help=f"Decoded files reconstructed per job (default: {FILES_PER_JOB})")
    parser.add_argument("--data-dir", type=str, default=DATA_DIR,
                        help=f"Directory of decoded .hipo input files (default: {DATA_DIR})")
    # AHDC Delta-phi K calibration: re-cook with an assumed beam position (mm) so
    # the modulation can be re-extracted (docs/ahdc_K_calibration_plan.md).
    parser.add_argument("--beam-x", type=float, default=0.0,
                        help="Assumed beam x (mm) for the AHDC KF beamline constraint (ALERT.beamX); 0 = production")
    parser.add_argument("--beam-y", type=float, default=0.0,
                        help="Assumed beam y (mm) for the AHDC KF beamline constraint (ALERT.beamY); 0 = production")
    parser.add_argument("--base-yaml", type=str, default=YAML_PATH,
                        help=f"Base recon YAML the ALERT beamX/beamY block is spliced into (default: {YAML_PATH})")
    parser.add_argument("--tag", type=str, default=None,
                        help="Suffix for the workflow/output name, to keep calibration passes separate "
                             "(default: auto 'bs_bx<X>_by<Y>' when a beam offset is given)")
    parser.add_argument("--submit", action="store_true",
                        help="After writing the JSON, submit it with 'swif2 import' + 'swif2 run' (run on ifarm)")
    args = parser.parse_args()

    data_files = read_data_files(args.data_dir)
    if not data_files:
        print("ERROR: No .hipo input files found. Exiting.")
        return

    n_events = -1
    if args.debug:
        data_files = data_files[:args.debug_count]
        n_events = args.debug_events
        print(f"DEBUG mode: {len(data_files)} file(s) total, {n_events} events/file, "
              f"partition={DEBUG_PARTITION}")

    if args.debug:
        resources = with_partition(recon_resources, DEBUG_PARTITION)
        resources["time_secs"] = DEBUG_TIME_SECS  # priority partition rejects the long production walltime
    else:
        resources = recon_resources

    batches = chunk(data_files, args.files_per_job)

    print(f"Generating 1 workflow: {len(batches)} job(s) "
          f"({args.files_per_job} files/job, {THREADS} threads)")

    out_dir = Path("json")
    out_dir.mkdir(parents=True, exist_ok=True)

    # K-calibration tagging: a beam offset (or explicit --tag) gets its own
    # workflow/output dir so passes never clobber each other on /volatile.
    use_beam = (args.beam_x != 0.0 or args.beam_y != 0.0 or args.tag is not None)
    tag = args.tag
    if tag is None and use_beam:
        tag = f"bs_bx{args.beam_x:g}_by{args.beam_y:g}".replace("-", "m").replace(".", "p")

    base_name: str = workflow_name_for()
    workflow_name = base_name
    if tag:
        workflow_name = f"{workflow_name}_{tag}"
    if args.debug:
        workflow_name = f"{workflow_name}_debug"
    # In debug mode keep outputs separate from production (…_debug/recon).
    output_recon: str = f"{OUTPUT_BASE}/{workflow_name}/recon"

    try:
        os.makedirs(output_recon, exist_ok=True)
    except OSError as exc:
        print(f"WARNING: could not create {output_recon} ({exc}). "
              f"Create it on the farm before submitting.")

    # Pick the recon YAML: production base, or a derived one carrying the assumed
    # beam position (ALERT beamX/beamY) for the AHDC Delta-phi K calibration.
    yaml_path = args.base_yaml
    if use_beam:
        yaml_path = f"{OUTPUT_BASE}/{workflow_name}/config.yaml"
        try:
            write_calibration_yaml(args.base_yaml, yaml_path, args.beam_x, args.beam_y)
        except OSError as exc:
            print(f"WARNING: could not write calibration YAML to {yaml_path} ({exc}). "
                  f"Generate it on the farm before submitting.")

    jobs: List[collections.OrderedDict] = []
    for i, batch in enumerate(batches):
        job_name: str = f"{base_name}_recon_{i}"
        # Each job gets its OWN run-clara output dir. run-clara writes a
        # per-run filelist.txt (and log/) into its -o directory; if all jobs
        # of a workflow share one recon dir on /volatile they clobber each
        # other's filelist.txt, and jobs then die with "Input files do not
        # exist" / "empty list of input files". A per-job subdir isolates it.
        job_output_recon: str = f"{output_recon}/job_{i}"
        command = build_command(batch, job_output_recon, yaml_path, n_events)
        jobs.append(create_recon_job(job_name, command, batch, len(batches), i, resources))

    workflow = collections.OrderedDict(
        {"name": workflow_name, "site": "jlab/enp", "max_dispatched": 1500, "jobs": jobs})

    output_path = out_dir / f"{workflow_name}.json"
    with output_path.open("w") as f:
        json.dump(workflow, f, **__JSONFORMAT)
    print(f"  {output_path} ({len(jobs)} jobs)")

    if args.submit:
        submit_workflow(workflow_name, output_path)

    print(f"Done: workflow JSON written under {out_dir}/")


if __name__ == "__main__":
    main()

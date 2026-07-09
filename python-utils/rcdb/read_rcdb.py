from rcdb import RCDBProvider

# Connect (MySQL or SQLite)
db = RCDBProvider("mysql://rcdb@clasdb-farm.jlab.org/rcdb")

# # Select run values with a query
# table = db.select_values(["beam_current", "event_count"],
#                          "@is_production and beam_current > 100",
#                          run_min=30000, run_max=31000)

table = db.select_values(["beam_current", "torus_scale"],
                         run_min=21000, run_max=23065)

# for row in table:
#     run_number, beam_current, event_count = row
#     print(f"{run_number},{beam_current},{event_count}")

with open("./output/results.csv", "w") as f :
    # write header
    f.write("run,beam,torus\n")
    for row in table:
        # write entry
        run_number, beam_current, event_count = row
        f.write(f"{run_number},{beam_current},{event_count}\n")
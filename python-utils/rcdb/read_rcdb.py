from rcdb import RCDBProvider

# Connect (MySQL or SQLite)
db = RCDBProvider("mysql://rcdb@clasdb-farm.jlab.org/rcdb")

#  # Example 1 with condition
# table = db.select_values(["beam_current", "torus_scale", "event_count"], search_str="event_count>1000",
#                          run_min=22222, run_max=23065)

# for row in table:
#     run_number, beam_current, torus_scale, event_count = row
#     print(f"{run_number},{beam_current},{torus_scale},{event_count}")


# Example 2 with no conditions
table = db.select_values(["beam_current", "torus_scale", "event_count"],
                         run_min=22222, run_max=23065)

with open("./output/results.csv", "w") as f :
    # write header
    f.write("run,beam,torus,nevents\n")
    for row in table:
        # write entry
        run_number, beam_current, torus, event_count = row
        f.write(f"{run_number},{beam_current},{torus},{event_count}\n")
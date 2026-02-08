import pandas as pd
from io import StringIO, BytesIO
import re
def parse_log(filename):
    with open(filename, "r") as f:
        lines = f.read().splitlines()

    # thermo data starts with a line containing "Step"
    # and ends just before a line containing "Loop"
    try:
        log_start_index = next(i for i, line in enumerate(lines) if line.strip().startswith("Step"))
        log_end_index = next(i for i, line in enumerate(lines) if "Loop" in line)
    except StopIteration:
        raise RuntimeError("Could not find Step/Loop markers in the log file")

    thermo_lines = [
        line
        for line in lines[log_start_index:log_end_index]
        if not line.strip().startswith("WARNING:")
    ]

    parsed_lines = [",".join(line.split()) for line in thermo_lines]

    csv_buffer = StringIO("\n".join(parsed_lines))

    return pd.read_csv(csv_buffer)


def parse_chunk(filename):
    with open(filename, "r") as f:
        lines = f.read().splitlines()

    # extract the timestamp from line 4 (index 3)
    end_timestamp = lines[3].split()[0]

    # extract headers from line 3 after the "# " (index 2)
    raw_headers = lines[2].split()[1:]

    # clean the headers: remove "c_", "v_", "[" and "]"
    cleaned_headers = [re.sub(r"c_|v_|\[|\]", "", h) for h in raw_headers]

    final_headers = ["timestep"] + cleaned_headers

    # parse data rows
    parsed_rows = []
    for line in lines[4:]:
        row = line.split()
        parsed_rows.append(",".join([end_timestamp, *row]))

    raw_csv = "\n".join([",".join(final_headers)] + parsed_rows)

    csv_buffer = StringIO(raw_csv)

    return pd.read_csv(csv_buffer)

def parse_dump(filename):
    rows = []

    with open(filename, "r") as f:
        while True:
            line = f.readline()
            if not line:
                break

            # ITEM: TIMESTEP
            if line.strip() == "ITEM: TIMESTEP":
                timestep = int(f.readline().strip())

                # ITEM: NUMBER OF ATOMS
                f.readline()
                num_atoms = int(f.readline().strip())

                # ITEM: BOX BOUNDS ...
                f.readline()
                f.readline()
                f.readline()
                f.readline()

                # ITEM: ATOMS ...
                atom_header = f.readline().split()[2:]
                if not rows:
                    column_names = ["timestep"] + atom_header

                for _ in range(num_atoms):
                    line = f.readline()
                    rows.append([timestep, *line.split()])

    return pd.DataFrame(rows, columns=column_names)

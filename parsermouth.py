import pandas as pd
from io import StringIO, BytesIO
import re
import warnings
from contextlib import contextmanager

@contextmanager
def open_input(source, mode='r'):
    """
    Handles file paths, StringIO, and BytesIO.
    Mode 'r' for text, 'rb' for bytes.
    """
    if isinstance(source, (str, bytes)) and not isinstance(source, (StringIO, BytesIO)):
        # If it's a string/bytes path, open the file
        with open(source, mode) as f:
            yield f
    else:
        # If it's already a stream (StringIO or BytesIO), just yield it
        yield source

def read_log_general(filename_or_fileobj, is_binary=False):
    mode = 'rb' if is_binary else 'r'

    with open_input(filename_or_fileobj, mode=mode) as f:
        raw_data = f.read()

    # If the data is bytes, decode it to a string first
    if isinstance(raw_data, bytes):
        content = raw_data.decode('utf-8')
    else:
        content = raw_data

    lines = content.splitlines()

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

    full_text = "\n".join(thermo_lines)
    return pd.read_csv(StringIO(full_text), sep=r'\s+')

def parse_log(filename):
    warnings.warn(
        "`parse_log(filename)` function is deprecated, use `read_log_general(filename_or_fileobj, is_binary=False)` instead.",
        DeprecationWarning,
        stacklevel=2
    )
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

def read_chunk_general(filename_or_fileobj, is_binary=False):
    mode = 'rb' if is_binary else 'r'

    with open_input(filename_or_fileobj, mode=mode) as f:
        raw_data = f.read()

    content = raw_data.decode('utf-8') if isinstance(raw_data, bytes) else raw_data
    lines = [l.strip() for l in content.splitlines() if l.strip() != ""]

    # ---- Extract headers ONCE ----
    header_line = None
    for l in lines:
        if l.startswith("# Chunk "):
            header_line = l
            break

    if header_line is None:
        raise ValueError("Could not find header line")

    raw_headers = header_line.split()[1:]
    cleaned_headers = [re.sub(r"c_|v_|\[|\]", "", h) for h in raw_headers]

    final_headers = ["timestep"] + cleaned_headers

    all_dfs = []
    i = 0
    n_lines = len(lines)

    while i < n_lines:

        line = lines[i]

        # skip comments
        if line.startswith("#"):
            i += 1
            continue

        parts = line.split()

        # ---- detect timestep line ----
        if len(parts) >= 2 and parts[0].isdigit():
            timestep = int(parts[0])
            n_chunks = int(parts[1])

            # ---- read chunk rows ----
            data_start = i + 1
            data_end = data_start + n_chunks

            data_body = "\n".join(lines[data_start:data_end])

            df = pd.read_csv(
                StringIO(data_body),
                sep=r"\s+",
                names=final_headers[1:]
            )

            df.insert(0, "timestep", timestep)

            all_dfs.append(df)

            i = data_end
        else:
            i += 1

    return pd.concat(all_dfs, ignore_index=True)

def parse_chunk(filename):
    warnings.warn(
        "`parse_chunk(filename)` function is deprecated, use `read_chunk_general(filename_or_fileobj, is_binary=False)` instead.",
        DeprecationWarning,
        stacklevel=2
    )
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

def read_dump_general(filename_or_fileobj, is_binary=False):
    mode = 'rb' if is_binary else 'r'

    with open_input(filename_or_fileobj, mode=mode) as f:
        raw_data = f.read()

    if isinstance(raw_data, bytes):
        content = raw_data.decode('utf-8')
    else:
        content = raw_data

    # setup iterator for processing
    lines = content.splitlines()
    line_iter = iter(lines)
    all_data = []
    headers = None

    # process blocks
    for line in line_iter:
        line = line.strip()

        if line == "ITEM: TIMESTEP":
            timestep = int(next(line_iter).strip())

            # skip ITEM: NUMBER OF ATOMS + read count
            next(line_iter)
            num_atoms = int(next(line_iter).strip())

            # skip ITEM: BOX BOUNDS (4 lines)
            for _ in range(4): next(line_iter)

            # parse ITEM: ATOMS header
            atom_header_line = next(line_iter)
            if headers is None:
                headers = ["timestep"] + atom_header_line.split()[2:]

            # read N atoms
            for _ in range(num_atoms):
                atom_data = next(line_iter).split()
                # convert to numeric immediately to prevent Object dtypes
                all_data.append([timestep] + [float(x) for x in atom_data])

    df = pd.DataFrame(all_data, columns=headers)

    # recast columns to integer
    potential_ints = ['timestep', 'id', 'type']
    cast_map = {col: int for col in potential_ints if col in df.columns}

    return df.astype(cast_map)

def parse_dump(filename):
    warnings.warn(
        "`parse_dump(filename)` function is deprecated, use `parse_dump(filename_or_fileobj, is_binary=False)` instead.",
        DeprecationWarning,
        stacklevel=2
    )
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

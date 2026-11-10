from io import StringIO, BytesIO, TextIOWrapper
from contextlib import contextmanager
from itertools import islice

@contextmanager
def open_input(source, mode='r'):
    """
    Handles file paths, StringIO, and Bytes
    Mode 'r' for text, 'rb' for bytes.
    """
    if isinstance(source, (str, bytes)) and not isinstance(source, (StringIO, BytesIO)):
        # If it's a string/bytes path, open the file
        with open(source, mode) as f:
            yield f
    else:
        # If it's already a stream (StringIO or BytesIO), just yield it
        yield source

def parse_chunk(src, dst, is_binary=False):
    """Convert a LAMMPS fix ave/chunk output file to CSV.

    src: file name or file object (given to open_input)
    dst: file object open for text writing
    """
    header_done = False

    with open_input(src, mode="rb" if is_binary else "r") as raw:
        f = TextIOWrapper(raw, encoding="utf-8") if is_binary else raw
        lines = (s for s in map(str.strip, f) if s)     # skip blank lines

        for line in lines:
            if line.startswith("#"):
                if not header_done and line.startswith("# Chunk "):
                    cols = line.split()[1:]             # drop the "#" token
                    dst.write("timestep," + ",".join(cols) + "\n")
                    header_done = True
                continue

            parts = line.split()
            if len(parts) >= 2 and parts[0].isdigit():
                if not header_done:
                    raise ValueError("Data found before header line")
                prefix = parts[0] + ","
                n_chunks = int(parts[1])

                rows = [prefix + ",".join(l.split()) for l in islice(lines, n_chunks)]
                if len(rows) != n_chunks:
                    raise ValueError(
                        f"Truncated block: expected {n_chunks} rows, got {len(rows)}"
                    )
                if rows:
                    dst.write("\n".join(rows) + "\n")

    if not header_done:
        raise ValueError("Could not find header line")

def parse_dump(src, dst, is_binary=False):
    """Convert a LAMMPS text dump to CSV.

    src: file name or file object (given to open_input)
    dst: file object open for text writing
    """
    header_done = False

    with open_input(src, mode="rb" if is_binary else "r") as raw:
        f = TextIOWrapper(raw, encoding="utf-8") if is_binary else raw

        for line in f:
            if not line.startswith("ITEM: TIMESTEP"):
                continue

            prefix = next(f).strip() + ","
            next(f)                       # ITEM: NUMBER OF ATOMS
            n_atoms = int(next(f))
            for _ in islice(f, 4):        # ITEM: BOX BOUNDS + 3 bound lines
                pass

            atom_header = next(f).split()[2:]   # drop "ITEM: ATOMS"
            if not header_done:
                dst.write("timestep," + ",".join(atom_header) + "\n")
                header_done = True

            rows = [prefix + ",".join(l.split()) for l in islice(f, n_atoms)]
            if len(rows) != n_atoms:
                raise ValueError(f"Truncated frame: expected {n_atoms} atoms, got {len(rows)}")
            dst.write("\n".join(rows) + "\n")

    if not header_done:
        raise ValueError("Could not find header line")

def parse_log(src, dst, is_binary=False):
    """Write the first thermo block of a LAMMPS log to CSV.

    src: file name or file object (given to open_input)
    dst: file object open for text writing
    """
    started = False

    with open_input(src, mode="rb" if is_binary else "r") as raw:
        f = TextIOWrapper(raw, encoding="utf-8") if is_binary else raw

        for line in f:
            if not started:
                if not line.lstrip().startswith("Step"):
                    continue
                started = True
            elif "Loop" in line:
                return                      # end of the thermo block
            elif line.lstrip().startswith("WARNING:"):
                continue

            dst.write(",".join(line.split()) + "\n")

    raise RuntimeError("Could not find Step/Loop markers in the log file")

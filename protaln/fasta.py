"""Minimal FASTA I/O."""


def parse_fasta(text):
    """Return a list of (header, sequence) tuples."""
    records, header, parts = [], None, []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(parts)))
            header, parts = line[1:].strip(), []
        elif header is not None:
            parts.append(line.replace(" ", "").upper())
        else:
            raise ValueError("FASTA data must start with a '>' header")
    if header is not None:
        records.append((header, "".join(parts)))
    return records


def read_fasta(path):
    with open(path) as handle:
        return parse_fasta(handle.read())


def write_fasta(records, handle, width=60):
    for header, seq in records:
        handle.write(">%s\n" % header)
        for i in range(0, len(seq), width):
            handle.write(seq[i:i + width] + "\n")

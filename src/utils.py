import pathlib as p


def concat_fastas(dst: p.Path, *fastas: p.Path):
    """Concat a list of fasta files into one fasta.

    Arguments
    ---------
    dst : Path
        Destination name for the resulting fasta.
    *fastas : Path
        Fasta files that should be concated together.
    """
    contents = [file.read_text() for file in fastas]
    with dst.open('w') as output:
        output.writelines('\n'.join(contents))

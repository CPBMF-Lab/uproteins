import pathlib as p


def concat_fastas(out: p.Path, *args: p.Path):
    contents = [file.read_text() for file in args]
    with out.open('w') as output:
        output.writelines('\n'.join(contents))

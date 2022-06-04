import os, sys

for file in sys.argv[1:]:
    with open("TEMPLATE.tex", "r") as f:
        source = f.read().replace("FILENAME", os.path.abspath(file))

        filename = file.split(".tex")[0].split("/")[-1]

        os.system(f"mkdir .{filename}")

        with open(f".{filename}/source.tex", "w") as s:
            s.write(source)

        os.system(
            f"pdflatex -aux-directory=.{filename} -output-directory=.{filename} .{filename}/source.tex "
        )
        os.system(f"pdf2svg .{filename}/source.pdf {filename}.svg")
        os.system(f"rm -rf .{filename}")

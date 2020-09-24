"""Script for parsing all stability plots to a markdown page."""
import os

DOC_ROOT = os.path.dirname(__file__)


def main():
    """Parses img dir and moves images to index.md."""
    text = "# Supporting figures and analysis code for our two-nucleon calculation with sLapH method\n\n"

    for root, _, files in os.walk("imgs"):
        section = " ".join(root.replace("imgs/", "").split("_")).capitalize()
        text += section + "\n\n"
        for ff in files:
            if not ff.endswith(".svg"):
                continue
            name = (
                " ".join(ff.replace(".svg", "").split("_"))
                .capitalize()
                .replace("Nn", "NN")
            )
            address = os.path.join(root, ff)
            text += f"## {name}\n![{name} plot]({address})\n\n"

    with open(os.path.join(DOC_ROOT, "index.md"), "w") as out:
        out.write(text)


if __name__ == "__main__":
    main()

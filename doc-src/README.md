# nn_c103_qcotd_swave_only


## About the webpage
The webpage is build with [MkDocs](https://www.mkdocs.org)

**The following commands are assumed to be run in the repo root dir**

Install:
```
pip install -r requirements.txt
```

Run locally
```
mkdocs serve
```

Update: Change the markdown files in `doc-src`.

Push:
```
mkdocs build
git add docs/
git commit -m "new version of page ..."
git push
```

## Updating the page

The repo root file `mkdocs.yml` specifies the webpage config.

The `doc-src` files (markdown) are the content of the pages.

## Converting pdfs to svgs:

See the `Makefile` for utilizing `inkscape`.

## Including svgs in the markdown

See the `index.md` file.

# Contributor guide

Contributions are welcome — bug reports, docs, and code. Please
[open an issue or pull request](https://github.com/hdrake/regionate/issues).

## Development environment

```bash
git clone https://github.com/hdrake/regionate.git
cd regionate
conda env create -f docs/environment.yml   # or: pip install -e .
conda activate docs_env_regionate
pip install -e .
```

Run the test suite with `pytest`.

## Example notebooks

Example notebooks live at the repository root in `examples/` so they run in
place, and are copied into `docs/examples/` at build time (see `DOC_NOTEBOOKS`
in `docs/conf.py`).

```{important}
Run notebooks locally and commit them **with their outputs**. The documentation
build does **not** execute notebooks (`nb_execution_mode = "off"` in
`docs/conf.py`), so whatever outputs you commit are exactly what readers see.
Re-run and re-commit a notebook whenever its code or results change.
```

## Building the docs

The documentation uses Sphinx + [Furo](https://pradyunsg.me/furo/) +
[myst-nb](https://myst-nb.readthedocs.io/). Build it locally exactly as CI and
Read the Docs do (warnings are errors):

```bash
pip install -r docs/requirements.txt
python -m sphinx -b html -W --keep-going docs docs/_build/html
```

Then open `docs/_build/html/index.html`.

# Configuration file for the Sphinx documentation builder.
# https://www.sphinx-doc.org/en/master/usage/configuration.html
#
# Shared "xeos stack" config for hdrake's x* packages: Sphinx + Furo + myst-nb,
# recursive-autosummary API reference, notebooks rendered from committed outputs
# (never executed at build time).

from importlib.metadata import version as get_version
from pathlib import Path
import shutil

# -- Project information -----------------------------------------------------
project = "regionate"                     # e.g. "regionate"
author = "Henri F. Drake"                       # e.g. "Henri F. Drake"
copyright = "2026, Henri F. Drake"                 # e.g. "2026, Henri F. Drake"

release = get_version(project)
version = ".".join(release.split(".")[:2])

master_doc = "index"

# -- General configuration ---------------------------------------------------
extensions = [
    # myst-nb bundles the MyST Markdown parser AND notebook support. Do NOT also
    # list myst_parser — a second registration errors.
    "myst_nb",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "sphinx_autodoc_typehints",
]

templates_path = ["_templates"]
# Do NOT exclude "generated" — that's the recursive-autosummary stub dir
# (docs/generated/, git-ignored) and Sphinx must read it into the doctree.
exclude_patterns = ["_build", "build", "Thumbs.db", ".DS_Store"]

# Example notebooks are authored with heading-level jumps (H1 -> H3); that trips
# myst-nb's structural linter. Suppress rather than edit committed notebooks.
suppress_warnings = ["myst.header"]

# -- Notebooks (myst-nb): render committed outputs, never execute on build ----
nb_execution_mode = "off"
myst_enable_extensions = ["dollarmath", "amsmath", "colon_fence"]

# -- API reference (recursive autosummary) -----------------------------------
# One page per object: the stub templates (_templates/autosummary/) drive
# :members: on the per-class/function pages, so we do NOT also set "members" as
# a global autodoc default (that double-documents and trips -W).
autosummary_generate = True
autosummary_imported_members = False
autodoc_default_options = {"show-inheritance": True}
autodoc_typehints = "description"
napoleon_numpy_docstring = True
napoleon_google_docstring = False
# Heavy import-time deps that need not be installed for autodoc to import the
# package. E.g. regionate: ["cartopy"]. Leave [] if the package pip-installs cleanly.
autodoc_mock_imports = []

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "myst-nb",
    ".ipynb": "myst-nb",
}

# -- Intersphinx: the package "family" web -----------------------------------
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "xarray": ("https://docs.xarray.dev/en/stable/", None),
    "xgcm": ("https://xgcm.readthedocs.io/en/stable/", None),
    "xeos": ("https://xeos.readthedocs.io/en/latest/", None),
    "sectionate": ("https://sectionate.readthedocs.io/en/latest/", None),
    "regionate": ("https://regionate.readthedocs.io/en/latest/", None),
    "xbudget": ("https://xbudget.readthedocs.io/en/latest/", None),
    "xwmt": ("https://xwmt.readthedocs.io/en/latest/", None),
    "xwmb": ("https://xwmb.readthedocs.io/en/latest/", None),
}

# -- Copy example notebooks from repo-root examples/ into docs/examples/ ------
HERE = Path(__file__).resolve()
DOCS_DIR = HERE.parent
REPO_ROOT = HERE.parents[1]
EXAMPLES_SRC = REPO_ROOT / "examples"
EXAMPLES_DST = DOCS_DIR / "examples"

# Only the notebooks referenced from the docs toctree (keeps the build light).
DOC_NOTEBOOKS = [
    "1_thickness_budget.ipynb",
    "2_advective_heat_convergence.ipynb",
    "3_Arctic_heat_CM4p25.ipynb",
    "4_bounded_by_named_sections.ipynb",
    "5_ECCO_LLC90_multiface_regions.ipynb",
    "6_idealized_corner_cases.ipynb",
]


def _sync_examples():
    if EXAMPLES_DST.exists():
        shutil.rmtree(EXAMPLES_DST)
    EXAMPLES_DST.mkdir(parents=True, exist_ok=True)
    for name in DOC_NOTEBOOKS:
        src = EXAMPLES_SRC / name
        if src.exists():
            shutil.copy2(src, EXAMPLES_DST / name)


_sync_examples()

# -- HTML output (Furo) ------------------------------------------------------
html_theme = "furo"
html_title = f"{project} {release}"
html_static_path = ["_static"]
html_css_files = ["custom.css"]  # styles the sidebar "View on GitHub" button

html_theme_options = {
    # `source_repository` is the repo URL the sidebar "View on GitHub" button
    # (docs/_templates/sidebar/brand.html) links to. Furo's own top-right
    # view/edit icons are disabled ("top_of_page_buttons": []): they point at the
    # page *source*, not the repo.
    "source_repository": "https://github.com/hdrake/regionate",    # e.g. "https://github.com/hdrake/regionate/"
    "source_branch": "main",
    "source_directory": "docs/",
    "top_of_page_buttons": [],
    # Per-package accent color (see README palette table).
    # NB: Furo prepends "--" to these keys itself — do NOT include it here, or
    # you get "----color-brand-primary" and the accent silently no-ops.
    "light_css_variables": {
        "color-brand-primary": "#b45309",
        "color-brand-content": "#b45309",
    },
    "dark_css_variables": {
        "color-brand-primary": "#f0a868",
        "color-brand-content": "#f0a868",
    },
    # A second GitHub link (repo home) in the footer.
    "footer_icons": [
        {
            "name": "GitHub",
            "url": "https://github.com/hdrake/regionate",
            "html": """
                <svg stroke="currentColor" fill="currentColor" stroke-width="0"
                     viewBox="0 0 16 16" width="1em" height="1em">
                  <path fill-rule="evenodd" d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0 0 16 8c0-4.42-3.58-8-8-8z"></path>
                </svg>
            """,
            "class": "",
        },
    ],
}

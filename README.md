# regionate

A package for creating xgcm-grid consistent regional masks and boundaries, leveraging its sibling package [`sectionate`](https://github.com/MOM6-community/sectionate).

[![PyPI](https://badge.fury.io/py/regionate.svg)](https://badge.fury.io/py/regionate)
[![Conda Version](https://img.shields.io/conda/vn/conda-forge/regionate)](https://anaconda.org/conda-forge/regionate)
[![Docs](https://readthedocs.org/projects/regionate/badge/?version=latest)](https://regionate.readthedocs.io/en/latest/)
[![License](https://img.shields.io/github/license/hdrake/regionate)](https://github.com/hdrake/regionate)

Quick Start Guide
-----------------

**For users: minimal installation within an existing environment**
```bash
conda install regionate
```

**For developers: installing from scratch using `conda`**
```bash
git clone git@github.com:hdrake/regionate.git
cd regionate
conda env create -f docs/environment.yml
conda activate docs_env_regionate
pip install -e .
python -m ipykernel install --user --name docs_env_regionate --display-name "docs_env_regionate"
jupyter-lab
```

# regionate

A package for creating xgcm-grid consistent regional masks and boundaries, leveraging its sibling package [`sectionate`](https://github.com/raphaeldussin/sectionate).

Quick Start Guide
-----------------

**For users: minimal installation within an existing environment**
```bash
pip install git+https://github.com/hdrake/regionate.git@main
```

**For developers: installing from scratch using `conda`**
```bash
git clone git@github.com:hdrake/regionate.git
cd regionate
conda env create -f ci/environment.yml
conda activate test_env_regionate
pip install -e .
python -m ipykernel install --user --name regionate --display-name "test_env_regionate"
jupyter-lab
```

# Grid-consistent region masks and boundary indices for ocean models

Leverages sibling package [`sectionate`](https://github.com/raphaeldussin/sectionate).

Quick Start Guide
-----------------

**Installing from scratch using `conda`**
```bash
git clone git@github.com:hdrake/regionate.git
cd regionate
conda env create -f ci/regionate.yml
conda activate regionate
pip install -e .
python -m ipykernel install --user --name regionate --display-name "regionate"
jupyter-lab
```

**Minimal installation within an existing environment**
```bash
pip install git+https://github.com/hdrake/regionate.git
```

# Installation

## Requirements

regionate is compatible with python 3 (>= version 3.11). It requires
[xgcm](https://xgcm.readthedocs.io/en/stable/) (>= 0.10.1, which ships the
bipolar north-fold boundary and the multi-tile `face_connections` padding fixes)
and the topology-driven API of its sibling package
[sectionate](https://github.com/MOM6-community/sectionate).

The topology-driven sectionate API is not on PyPI yet, so install it from its
development branch first:

```bash
pip install "sectionate @ git+https://github.com/hdrake/sectionate.git@topology-driven-neighbors"
```

## From conda-forge

```bash
conda install -c conda-forge regionate
```

## From PyPI

```bash
pip install regionate
```

This installs the latest release from [PyPI](https://pypi.python.org/pypi).

## From GitHub (development version)

regionate is under active development. To obtain the latest development version,
clone the [source repository](https://github.com/hdrake/regionate) and install
it with pip:

```bash
pip install git+https://github.com/hdrake/regionate.git
```

See the [Contributor Guide](contributing.md) for a full development setup. Users
are encouraged to [fork](https://help.github.com/articles/fork-a-repo/) regionate
and submit [issues](https://github.com/hdrake/regionate/issues) and
[pull requests](https://github.com/hdrake/regionate/pulls).

# How we publish to PyPI

Do this from the repo root each time we ship a new version.

## 1. Bump the version

Put the same new version in these three files. PyPI refuses a version that is already up there, so always move the number forward (`2.0.1`, `2.1.0`, whatever is next).

- `pyproject.toml`
- `setup.py`
- `pygamlab/__init__.py` (`__version__`)

## 2. Build

First time on a machine:

```bash
python3 -m pip install build twine
```

Then, from the repo root:

```bash
rm -rf dist build *.egg-info
python3 -m build
twine check dist/*
```

`twine check` should say both files pass. That creates `dist/` with a `.tar.gz` and a `.whl`. Do not commit `dist/`, `build/`, or `*.egg-info/`.

## 3. Upload

PyPI wants an API token, not the account password. Make one at https://pypi.org/manage/account/token/ for this project.

```bash
twine upload dist/*
```

When it asks:

- username: `__token__`
- password: the token (`pypi-...`)

When the upload finishes, the new version is on https://pypi.org/project/pygamlab/

Install it and check the version:

```bash
python3 -m pip install --upgrade pygamlab
python3 -c "import pygamlab; print(pygamlab.__version__)"
```

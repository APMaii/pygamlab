# How we build the docs

Do this from the repo root. Edit pages in `doc_preparing/`. The site that GitHub Pages publishes stays in `doc/`.

## 1. Save the current site

Copy the live `doc/` folder into `archives/` before you replace it. Use the version you are putting away, and pick a new folder name if that one already exists.

```bash
cp -a doc archives/doc_2.0.0
```

## 2. Build

First time on a machine:

```bash
python3 -m pip install sphinx furo nbsphinx sphinx-autodoc-typehints
python3 -m pip install -e .
```

Then, from the repo root:

```bash
rm -rf doc_build
sphinx-build -b html doc_preparing doc_build
```

Open `doc_build/index.html` and look through it. If the build fails, fix the page in `doc_preparing/` and run the same command again. `doc/` is still the old site at this point.

## 3. Put the new site in doc/

```bash
rm -rf doc
cp -a doc_build doc
cp doc/index.html doc/documentation.html
cp -a doc_preparing/extra/. doc/
rm -rf doc_build
```

`cp doc/index.html doc/documentation.html` keeps the API/docs homepage at `documentation.html`.

`doc_preparing/extra/` puts back the custom pages Sphinx does not make or should not replace: `index.html`, `index_en.html`, `index_fa.html`, and `quickstart_tutorial_fa.html`.

Do not commit `doc_build/`. Commit `doc/` and the new folder under `archives/`. A push to `main` publishes `doc/` through the GitHub Pages workflow.

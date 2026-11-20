#!/usr/bin/env bash
# Regenerate docs/tutorial.rst and docs/tutorial_files/ from the tutorial
# notebook. The notebook is executed first, so the page always shows the
# output of the installed pystackreg.
#
# Requires: jupyter nbconvert, ipykernel, matplotlib, scikit-image, tifffile
# and an installed pystackreg. Run from the docs/ directory.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

jupyter nbconvert --to notebook --execute --inplace \
    --ExecutePreprocessor.timeout=600 ../examples/notebooks/tutorial.ipynb
jupyter nbconvert ../examples/notebooks/tutorial.ipynb --to rst
mv ../examples/notebooks/tutorial.rst .
# "ipython3" is not a Pygments lexer unless IPython's Sphinx extension is loaded
sed -i 's/^\.\. code:: ipython3/.. code:: python/' tutorial.rst
rm -rf tutorial_files/
mv ../examples/notebooks/tutorial_files .

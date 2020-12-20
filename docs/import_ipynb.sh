jupyter nbconvert ../examples/notebooks/tutorial.ipynb --to rst
mv ../examples/notebooks/tutorial.rst .
rm -rf tutorial_files/
mv ../examples/notebooks/tutorial_files .

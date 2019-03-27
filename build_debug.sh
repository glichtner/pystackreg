rm -rf build
rm pystackreg/*.so
python setup_debug.py build_ext --force --inplace --debug

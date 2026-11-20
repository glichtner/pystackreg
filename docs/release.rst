Release procedure
=================

Releases are made by the maintainer. Nothing is published automatically from
ordinary commits: wheels are uploaded to PyPI only by the
``wheels-deploy.yml`` workflow when a tag ``v*`` is pushed.

1. **Prepare**

   * Make sure the ``ci.yml`` workflow is green on ``master``.
   * Update ``CHANGELOG.rst``: replace "(unreleased)" by the release date.
   * Set the version in ``pystackreg/version.py`` (the single source of the
     version number), e.g. ``0.3.0.dev0`` -> ``0.3.0``. With bump2version:
     ``bump2version release`` (see ``.bumpversion.cfg``).
   * If the numerical code changed, run the exact comparison on Linux
     x86-64::

         PYSTACKREG_EXACT_BASELINE=1 pytest tests/test_regression.py

2. **Build and check locally** in a clean clone::

       git clone https://github.com/glichtner/pystackreg && cd pystackreg
       python -m pip install build twine
       python -m build                      # sdist and wheel in dist/
       python -m twine check --strict dist/*
       tar -tzf dist/pystackreg-*.tar.gz | grep -E "LICENSE|NOTICE"
       unzip -l dist/pystackreg-*.whl | grep -E "LICENSE|NOTICE|turboreg"

   Install the wheel into a fresh environment and run the tests against it
   from outside the source tree::

       python -m venv /tmp/prs && /tmp/prs/bin/pip install dist/pystackreg-*.whl pytest tifffile scipy
       cd /tmp && PYSTACKREG_TEST_INSTALLED=1 /tmp/prs/bin/pytest /path/to/pystackreg/tests

3. **Tag**::

       git commit -am "Release 0.3.0"
       git tag -a v0.3.0 -m "pystackreg 0.3.0"
       git push origin master v0.3.0

   Pushing the tag starts ``wheels-deploy.yml``, which builds the wheels with
   cibuildwheel for all platforms, tests each wheel, builds the sdist, checks
   the metadata and uploads everything to PyPI.

4. **After the upload**

   * Create a GitHub release for the tag with the changelog section as
     description. If the repository is connected to Zenodo, this archives the
     release and assigns a DOI; add the DOI to ``CITATION.cff``.
   * conda-forge: the feedstock bot opens a pull request for the new version.
     For 0.3.0 the recipe must be edited by hand: ``license: GPL-3.0-or-later``,
     ``license_file`` (``LICENSE`` and ``NOTICE``), ``python >=3.10``.
   * Check that https://pystackreg.readthedocs.io built the new version.
   * Start the next development cycle: set the version to the next
     ``.dev0`` and add an "(unreleased)" section to the changelog.

Manual builds
-------------

``wheels-deploy.yml`` can also be started by hand (``workflow_dispatch``). It
then builds and tests all wheels and stores them as workflow artifacts
without uploading anything.

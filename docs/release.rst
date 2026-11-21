Release procedure
=================

Releases are made by the maintainer. Nothing is published automatically from
ordinary commits: wheels are uploaded to PyPI only by the
``wheels-deploy.yml`` workflow when a tag ``v*`` is pushed.

One-time setup
--------------

The upload uses PyPI's `trusted publishing
<https://docs.pypi.org/trusted-publishers/>`_: the workflow identifies itself
to PyPI with a short-lived token, and no API token is stored in the
repository. PyPI has to know the workflow once. On https://pypi.org, open
*Manage project* -> *Publishing* for pystackreg and add a GitHub publisher
with

* owner ``glichtner``, repository ``pystackreg``,
* workflow ``wheels-deploy.yml``,
* environment ``pypi``.

GitHub creates the environment ``pypi`` when the upload job runs for the first
time. A required reviewer can be added to it in the repository settings
(*Settings* -> *Environments*); the upload then waits for a confirmation.

Steps
-----

1. **Prepare**

   * Make sure the ``ci.yml`` workflow is green on ``master``.
   * Update ``CHANGELOG.rst``: replace "(unreleased)" by the release date and
     commit this change. bump2version (``pip install bump2version``) refuses
     to run with uncommitted changes.
   * Set the version in ``pystackreg/version.py`` (the single source of the
     version number) with ``bump2version release``, e.g. ``0.3.0.dev0`` ->
     ``0.3.0`` (see ``.bumpversion.cfg``). This creates a commit but no tag.
   * If the numerical code changed, run the exact comparison on Linux
     x86-64::

         PYSTACKREG_EXACT_BASELINE=1 pytest tests/test_regression.py

     All cases are bit-identical only on a machine like the one that
     generated the baseline; elsewhere the cases with synthetic images can
     differ in the last bits (see the validation page).

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

       git tag -a v0.3.0 -m "pystackreg 0.3.0"
       git push origin master v0.3.0

   Pushing the tag starts ``wheels-deploy.yml``, which builds the wheels with
   cibuildwheel for all platforms, tests each wheel, builds the sdist, checks
   the metadata and uploads everything to PyPI.

4. **After the upload**

   * Create a GitHub release for the tag with the changelog section, in
     Markdown, as description.
   * Archive the release on Zenodo. The repository is not connected to
     Zenodo, so this is done by hand: open the most recent version through
     the concept DOI https://doi.org/10.5281/zenodo.24425557, choose *New
     version*, replace the file by the source archive of the new tag
     (``pystackreg-X.Y.Z.zip`` from the GitHub release), update the version
     and the publication date, and publish. The publisher stays "Zenodo".
     All versions share the concept DOI, which ``CITATION.cff`` and the
     README cite, so nothing has to be changed in the repository.
   * conda-forge: the feedstock bot opens a pull request for the new version.
     For 0.3.0 the recipe must be edited by hand: ``license: GPL-3.0-or-later``,
     ``license_file`` (``LICENSE`` and ``NOTICE``), ``python >=3.10`` and the
     build requirements ``numpy >=2.0`` and ``setuptools >=77``.
   * Check that https://pystackreg.readthedocs.io built the new version.
   * Start the next development cycle: ``bump2version patch`` sets the
     version to the next ``.dev0`` (``0.3.0`` -> ``0.3.1.dev0``); add an
     "(unreleased)" section to the changelog.

Manual builds
-------------

``wheels-deploy.yml`` can also be started by hand (``workflow_dispatch``). It
then builds and tests all wheels and stores them as workflow artifacts
without uploading anything.

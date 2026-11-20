"""Shared fixtures and environment checks for the pystackreg test suite."""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent))

import pystackreg

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = Path(__file__).resolve().parent / "data"


def pytest_report_header(config):
    return [
        f"pystackreg {pystackreg.__version__} imported from "
        f"{Path(pystackreg.__file__).parent}",
        f"numpy {np.__version__}",
    ]


def pytest_sessionstart(session):
    """Refuse to test the source tree when an installed build is expected.

    CI sets ``PYSTACKREG_TEST_INSTALLED=1`` when it tests a wheel. The check
    guards against silently importing the package directory of the checkout
    (and a stale extension module that may lie in it) instead.
    """
    if os.environ.get("PYSTACKREG_TEST_INSTALLED") != "1":
        return
    package_dir = Path(pystackreg.__file__).resolve().parent
    if package_dir == REPO_ROOT / "pystackreg":
        raise pytest.UsageError(
            "PYSTACKREG_TEST_INSTALLED=1, but pystackreg was imported from the "
            f"source tree ({package_dir}) instead of an installed distribution"
        )


@pytest.fixture(scope="session")
def pc12():
    """Unregistered example stack (5 frames, 201 x 199, uint16)."""
    tifffile = pytest.importorskip("tifffile")
    return tifffile.imread(DATA_DIR / "pc12-unreg.tif")


@pytest.fixture(scope="session")
def baseline():
    """Results of pystackreg 0.2.8 (see ``generate_baseline.py``)."""
    with np.load(DATA_DIR / "baseline_v0.2.8.npz") as data:
        return {key: data[key] for key in data.files}

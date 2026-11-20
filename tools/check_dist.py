"""Check the contents and metadata of built distributions.

    python tools/check_dist.py dist/

Verifies for every wheel and sdist in the directory that the license files
are included, that the metadata declares the expected license, that wheels
contain the compiled extension and no test or build files, and that the sdist
contains everything needed to build and test the package.
"""

from __future__ import annotations

import sys
import tarfile
import zipfile
from email.parser import Parser
from pathlib import Path

LICENSE_EXPRESSION = "GPL-3.0-or-later"
LICENSE_FILES = ("LICENSE", "NOTICE")
SDIST_REQUIRED = (
    "pyproject.toml",
    "setup.py",
    "README.rst",
    "CHANGELOG.rst",
    "LICENSE",
    "NOTICE",
    "inc/TurboReg.h",
    "inc/matrix.h",
    "src/pymain.cpp",
    "src/TurboRegTransform.cpp",
    "pystackreg/__init__.py",
    "pystackreg/py.typed",
    "tests/test_regression.py",
    "tests/data/baseline_v0.2.8.npz",
    "tests/data/pc12-unreg.tif",
)
WHEEL_REQUIRED = (
    "pystackreg/__init__.py",
    "pystackreg/pystackreg.py",
    "pystackreg/_internal.py",
    "pystackreg/util/__init__.py",
    "pystackreg/version.py",
    "pystackreg/py.typed",
)
WHEEL_FORBIDDEN = ("tests/", "conftest.py", "inc/", "src/", "docs/")


def check_metadata(text: str, errors: list[str]) -> None:
    meta = Parser().parsestr(text)
    if meta.get("License-Expression") != LICENSE_EXPRESSION:
        errors.append(f"License-Expression is {meta.get('License-Expression')!r}")
    declared = meta.get_all("License-File") or []
    for name in LICENSE_FILES:
        if name not in declared:
            errors.append(f"License-File {name} not declared in the metadata")
    if any(c.startswith("License ::") for c in meta.get_all("Classifier") or []):
        errors.append("legacy license classifier present")
    if not meta.get("Requires-Python"):
        errors.append("Requires-Python is missing")


def check_wheel(path: Path) -> list[str]:
    errors: list[str] = []
    with zipfile.ZipFile(path) as zf:
        names = zf.namelist()
        for required in WHEEL_REQUIRED:
            if required not in names:
                errors.append(f"missing {required}")
        if not any(
            n.startswith("pystackreg/turboreg.") and n.endswith((".so", ".pyd"))
            for n in names
        ):
            errors.append("compiled extension pystackreg/turboreg.* is missing")
        for name in names:
            if name.startswith(WHEEL_FORBIDDEN) or "/conftest.py" in name:
                errors.append(f"unexpected file {name}")
        dist_info = next(n.split("/")[0] for n in names if ".dist-info/" in n)
        for name in LICENSE_FILES:
            if f"{dist_info}/licenses/{name}" not in names:
                errors.append(f"{name} is not in {dist_info}/licenses/")
        check_metadata(zf.read(f"{dist_info}/METADATA").decode(), errors)
    return errors


def check_sdist(path: Path) -> list[str]:
    errors: list[str] = []
    with tarfile.open(path) as tf:
        names = tf.getnames()
        root = names[0].split("/")[0]
        members = {n[len(root) + 1 :] for n in names}
        for required in SDIST_REQUIRED:
            if required not in members:
                errors.append(f"missing {required}")
        for name in members:
            if name.endswith((".so", ".pyd", ".pyc", ".o")):
                errors.append(f"binary file {name}")
        pkg_info = tf.extractfile(f"{root}/PKG-INFO")
        assert pkg_info is not None
        check_metadata(pkg_info.read().decode(), errors)
    return errors


def main(directory: str) -> int:
    files = sorted(Path(directory).glob("*.whl")) + sorted(
        Path(directory).glob("*.tar.gz")
    )
    if not files:
        print(f"no distributions found in {directory}")
        return 1
    failed = False
    for path in files:
        errors = check_wheel(path) if path.suffix == ".whl" else check_sdist(path)
        print(f"{'FAIL' if errors else 'ok  '} {path.name}")
        for error in errors:
            print(f"       {error}")
        failed = failed or bool(errors)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "dist"))

# -*- coding: utf-8 -*-

import os
from setuptools import setup, Extension, find_packages

# cannot directly import because __init__.py imports pystackreg which imports the
# compiled plugin, which is not available before setup.py is run
__version__ = ""  # placeholder for linters
exec(open("pystackreg/version.py").read())


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


class get_numpy_include(object):
    def __str__(self):
        import numpy

        return numpy.get_include()


if __name__ == "__main__":
    readme = read("README.rst")
    change = read("CHANGELOG.rst")

    setup(
        name="pystackreg",
        description="Image registration tool (python implementation of the ImageJ/FIJI "
        "Plugin TurboReg/StackReg)",
        long_description="\n\n".join([readme, change]),
        version=__version__,
        author="Gregor Lichtner (python/C++ port);"
        "TurboReg Author: Philippe Thévenaz, Biomedical Imaging Group,"
        "Swiss Federal Institute of Technology Lausanne",
        url="https://github.com/glichtner/pystackreg",
        packages=find_packages("."),
        ext_modules=[
            Extension(
                "pystackreg.turboreg",
                [
                    "src/pymain.cpp",
                    "src/TurboReg.cpp",
                    "src/TurboRegMask.cpp",
                    "src/TurboRegImage.cpp",
                    "src/TurboRegTransform.cpp",
                    "src/TurboRegPointHandler.cpp",
                ],
                extra_compile_args=["-std=c++11"],
                include_dirs=["inc/", get_numpy_include()],
            )
        ],
        setup_requires=["numpy"],
        install_requires=["numpy", "tqdm"],
        classifiers=[
            # complete list: http://pypi.python.org/pypi?%3Aaction=list_classifiers
            "Development Status :: 5 - Production/Stable",
            "Intended Audience :: Science/Research",
            "License :: Free To Use But Restricted",
            "Operating System :: Microsoft :: Windows",
            "Operating System :: Unix",
            "Operating System :: POSIX",
            "Programming Language :: C++",
            "Programming Language :: Python",
            "Programming Language :: Python :: 2.7",
            "Programming Language :: Python :: 3.5",
            "Programming Language :: Python :: 3.6",
            "Programming Language :: Python :: 3.7",
            "Programming Language :: Python :: 3.8",
            "Programming Language :: Python :: 3.9",
            "Topic :: Scientific/Engineering :: Image Processing",
            "Topic :: Utilities",
        ],
    )

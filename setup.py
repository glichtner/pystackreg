# -*- coding: utf-8 -*-

import setuptools
from distutils.core import setup, Extension
import numpy.distutils.misc_util
import os


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


exec(open("pystackreg/version.py").read())


setup(
    name="pystackreg",
    description="Python implementation of the ImageJ/FIJI Plugin TurboReg/StackReg",
    long_description=read("README.rst"),
    version=__version__,
    author="Gregor Lichtner (python/C++ port); TurboReg Author: Philippe Thévenaz, Biomedical Imaging Group, Swiss Federal Institute of Technology Lausanne",
    url="https://github.com/glichtner/pystackreg",
    packages=["pystackreg"],
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
        )
    ],
    include_dirs=["inc/"] + numpy.distutils.misc_util.get_numpy_include_dirs(),
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

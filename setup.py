# -*- coding: utf-8 -*-

import setuptools
from distutils.core import setup, Extension
import numpy.distutils.misc_util
import os

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

exec(open('pystackreg/version.py').read())

setup(
	name="pystackreg",
	description='Python implementation of the ImageJ/FIJI Plugin TurboReg/StackReg',
	long_description=read('README.rst'),
	version=__version__,
	author='Gregor Lichtner (python/C++ port); TurboReg Author: Philippe Thévenaz, Biomedical Imaging Group, Swiss Federal Institute of Technology Lausanne',
	url='https://bitbucket.org/glichtner/pystackreg',
	#package_dir={'pystackreg':'pysrc'},
	packages=['pystackreg'],
	ext_modules=[
		Extension("pystackreg.turboreg", [
			"src/pymain.cpp",
			#"src/PyStackReg.cpp",
			"src/TurboReg.cpp",
			"src/TurboRegMask.cpp",
			"src/TurboRegImage.cpp",
			"src/TurboRegTransform.cpp",
			"src/TurboRegPointHandler.cpp"],
			)],
	include_dirs=['inc/'] + numpy.distutils.misc_util.get_numpy_include_dirs(),
	install_requires = ['numpy', 'tqdm']
) 


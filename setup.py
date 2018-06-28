# -*- coding: utf-8 -*-

from distutils.core import setup, Extension
import numpy.distutils.misc_util

setup(
	name="pystackreg",
	description='Python implementation of the ImageJ/FIJI Plugin TurboReg/StackReg',
	version="0.1",
	author='Gregor Lichtner (python/C++ port); TurboReg Author: Philippe Thévenaz, Biomedical Imaging Group, Swiss Federal Institute of Technology Lausanne',
	package_dir={'pystackreg' : 'pysrc'},
	packages=['pystackreg'],
    ext_modules=[
		Extension("pystackreg.stackreg", [
			"src/pymain.cpp",
			#"src/PyStackReg.cpp",
			"src/TurboReg.cpp",
			"src/TurboRegMask.cpp",
			"src/TurboRegImage.cpp",
			"src/TurboRegTransform.cpp",
			"src/TurboRegPointHandler.cpp"])],
    include_dirs=['inc/'] + numpy.distutils.misc_util.get_numpy_include_dirs(),
) 


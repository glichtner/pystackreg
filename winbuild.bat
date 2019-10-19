CALL c:\ProgramData\Anaconda2\scripts\activate py27
rd /s/q build
python setup.py build
python setup.py bdist_wheel

CALL c:\ProgramData\Anaconda2\scripts\activate py36
rd /s/q build
python setup.py build
python setup.py bdist_wheel

CALL c:\ProgramData\Anaconda2\scripts\activate py37
rd /s/q build
python setup.py build
python setup.py bdist_wheel

set CONDA_FORCE_32BIT=1

CALL c:\ProgramData\Anaconda2\scripts\activate py27_32bit
rd /s/q build
python setup.py build
python setup.py bdist_wheel

CALL c:\ProgramData\Anaconda2\scripts\activate py36_32bit
rd /s/q build
python setup.py build
python setup.py bdist_wheel

CALL c:\ProgramData\Anaconda2\scripts\activate py37_32bit
rd /s/q build
python setup.py build
python setup.py bdist_wheel

set CONDA_FORCE_32BIT=
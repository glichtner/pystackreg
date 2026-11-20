Authors
=======

Original method and software
----------------------------

* Registration method: Philippe Thévenaz, Urs E. Ruttimann and Michael Unser,
  *A Pyramid Approach to Subpixel Registration Based on Intensity*, IEEE
  Transactions on Image Processing 7(1), 27-41 (1998).
* Developer of the original ImageJ plugins TurboReg/StackReg: Philippe Thévenaz,
  Biomedical Imaging Group, EPFL (https://bigwww.epfl.ch/thevenaz/)

pystackreg
----------

* Developer and maintainer of this package (port of the Java code from
  P. Thévenaz to C++, Python interface and all further functionality):
  Gregor Lichtner `@glichtner <https://github.com/glichtner>`_

Contributors
............

* Sebastian Oltmanns - progress callback
  (`PR #1 <https://github.com/glichtner/pystackreg/pull/1>`_)
* Maksim Rakitin - license file in the source distribution
  (`PR #4 <https://github.com/glichtner/pystackreg/pull/4>`_)
* Florian Rey - installation without a preinstalled NumPy
  (`PR #10 <https://github.com/glichtner/pystackreg/pull/10>`_)
* Oren Amsalem - ``tqdm.auto`` progress bars, parallel stack registration with
  a process pool, documentation
  (`PR #11 <https://github.com/glichtner/pystackreg/pull/11>`_,
  `PR #16 <https://github.com/glichtner/pystackreg/pull/16>`_,
  `PR #17 <https://github.com/glichtner/pystackreg/pull/17>`_)
* Patrick Schuenke - fix for the ``axis`` argument with ``reference="mean"``
  (`PR #26 <https://github.com/glichtner/pystackreg/pull/26>`_)
* `@camilo1729 <https://github.com/camilo1729>`_ - proposal to release the GIL
  in the native code
  (`issue #14 <https://github.com/glichtner/pystackreg/issues/14>`_)

Thanks to everyone who reported problems and suggested improvements in the
`issue tracker <https://github.com/glichtner/pystackreg/issues>`_.

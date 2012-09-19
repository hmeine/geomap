This repository contains ten years worth of GeoMap code.

History
-------

The code was ported from boost::python to boost::python V2, from
vigra/interactive to vigra/interactive2, partially to
vigra/interactive3, and finally (at least much of it) to the current,
released vigranumpy bindings.  Furthermore, the GUI stuff was ported
from Qt3 to Qt4.

Status
------

However, the latest porting effort to vigranumpy (and maybe even to
Qt4) is not finished, i.e. not all code may be ported or tested yet.

Furthermore, during the preparation of this Git repository to be
published, a lot of cruft has been removed, and many subdirectories
were introduced. (Imagine all files in one directory, intertwined with
at least the same number of files from one-off experiments and student
thesises, and you are close to the status I was looking at this
morning.)

Compilation and running will *fail* at the moment, because the build
system, includes and imports have not been adapted yet.  Nevertheless,
we have a public repository now, which at last makes collaboration
possible! :-)

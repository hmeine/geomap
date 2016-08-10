This repository contains ten years worth of GeoMap code.

History
-------

The code was ported from boost::python to boost::python V2, from
vigra/interactive to vigra/interactive2, partially to
vigra/interactive3, and finally (at least much of it) to the current,
released vigranumpy bindings.  Furthermore, the GUI stuff was ported
from Qt3 to Qt4.  (And we could soon start the next port to Qt5.)

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

Compilation was fixed, but the python files have not been updated and tested yet.
Nevertheless, we have a public repository now, which at last makes
collaboration possible! :-)

Structure
---------

Interesting frontend bits that come to my mind (many of the often
needed and more or less stable functions and classes are even
documented):

* general helper functions are in maputils.py
* sub-pixel watersheds are initialized from within maputils.py
* crack-edge maps (and mid-cracks etc.) are initialized using crackConvert.py
* display GeoMaps on top of images using mapdisplay.py
* (Constrained) Delaunay Triangulation is in delaunay.py
* Chordal Axis Transform (CAT) is in delaunay.py, too
* a simple segmentation GUI / frontend is in workplace.py


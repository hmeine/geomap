The GeoMap is a unified representation for segmentation results, which
comprises both their geometry and topology.  The GeoMap is intended to be used
not only for the *final* segmentation result, but in particular for intermediate
working states of a complex segmentation process chain.  Hence, it offers
operations that modify the tesselation, for instance merging two neighboring
regions.  This particular implementation is using a polygonal representation,
and consists of a core C++ library, Python bindings, and more complex operations
(and GUI components) being implemented in Python only.

A list of corresponding scientific publications can be found at:
http://www.citeulike.org/user/hans_meine/tag/geomap

This repository contains ten years worth of GeoMap code.

Code History
------------

The code was ported from boost::python to boost::python V2, from
vigra/interactive to vigra/interactive2, partially to
vigra/interactive3, and finally (at least much of it) to the current,
released vigranumpy bindings.  Furthermore, the GUI stuff was ported
from Qt3 to Qt4.  (And we could soon start the next port to Qt5.)
The code was initially managed with CVS, then in Subversion, Mercurial,
and now Git.  It was always tried to preserve the history during
conversion, but of course chances are that it suffered.

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

License
-------
Licensed under `MIT license <LICENSE.txt>`_.

The current CMakeLists are based on the VIGRA infrastructure.  Simply
place this stuff under [VIGRA]/vigranumpy/private (directly or in a
subdirectory referenced by an appropriate CMakeLists.txt) and it will
be compilable via VIGRA's build system.

ATM, you seem to need "make geomap" since it is not part of the
default build.

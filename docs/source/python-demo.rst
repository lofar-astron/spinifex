.. _python-demo-label:
Python demo
===========
Example script how to get ionospheric rotation measures for a give source,
observer location and array of times. This script uses the default ionospheric settings,
which are good for most purposes.

.. literalinclude:: examples/example_get_rm_from_skycoord.py

You can also, for a given time, map the ionospheric rotation measure sky above the observer
by giving a set of azimuth and elevation coordinates:

.. literalinclude:: examples/example_get_rm_from_altaz.py

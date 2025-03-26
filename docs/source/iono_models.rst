==================
Ionospheric Models
==================


Spinifex has implemented different ionospheric models:

    * :ref:`ionex`
    * :ref:`ionex_iri`
    * :ref:`tomion`

.. _ionex:

ionex
---------------------
The fastest method is to use IONEX datafiles that are freely available from various online servers
and deliver global ionospheric models (GIM) with limited accuracy. A single layer ionospheric model is assumed at a
user specified height.

The ionex model comes with the following options:

.. _ionex_iri:

ionex_iri
---------------------
A more advanced method uses the integrated total electron content (TEC) from the IONEX files, but also includes
a density profile at user specified heights. The most important advantage of using the density profile
is a better estimate of the plasmaspheric contribution to the TEC. This avoids to a large extent the observed
overestimation of ionospheric Faraday rotation when a single layer is assumed.

The ionex_iri model comes with the following options:

.. _tomion:

tomion
---------------------
The 2 layer ionospheric model provided by UPC-IonSat. To be implemented soon...

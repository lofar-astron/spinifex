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

* height :  [u.Quantity] altitude of the single layer ionosphere. Default is 350 km
* server: [str] Server to download from, by default "cddis". Must be a supported server.
* prefix : [str] Analysis centre prefix, by default "cod". Must be a supported analysis centre.
* time_resolution : [u.Quantity] Time resolution, by default None, will default to the server time resolution.
* solution : Solution type, by default "final", must be "final" or "rapid".
* output_directory : [Path] Output directory path, by default None, will default to "ionex_files" in the current working directory.


.. _ionex_iri:

ionex_iri
---------------------
A more advanced method uses the integrated total electron content (TEC) from the IONEX files, but also includes
a normalised electron density profile from the international reference ionosphere (IRI). The most important advantage
of using the density profile
is a better estimate of the plasmaspheric contribution to the TEC. This avoids to a large extent the observed
overestimation of ionospheric Faraday rotation when a single layer is assumed.

The ionex_iri model comes with the following options:

* height :  [u.Quantity] altitude of the single layer ionosphere. Default is 350 km
* server: [str] Server to download from, by default "cddis". Must be a supported server.
* prefix : [str] Analysis centre prefix, by default "cod". Must be a supported analysis centre.
* time_resolution : [u.Quantity] Time resolution, by default None, will default to the server time resolution.
* solution : Solution type, by default "final", must be "final" or "rapid".
* output_directory : [Path] Output directory path, by default None, will default to "ionex_files" in the current working directory.



.. _tomion:

tomion
---------------------
The 2 layer ionospheric model provided by UPC-IonSat. To be implemented soon...

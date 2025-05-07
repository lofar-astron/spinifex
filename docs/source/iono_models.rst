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
* correct_uqrg_rms : [bool] The rms values for UQR maps are overestimated (see `Zhao et al. <https://link.springer.com/article/10.1007/s00190-021-01487-8>`_), use this flag to correct these. default=True

At the moment three hosts of ionex data are supported:

* NASA CDDIS - `cddis <https://cddis.nasa.gov/archive/gnss/products/ionex>`_,
* Barcelona Tech UPC Chapman - `chapman <http://chapman.upc.es/tomion/rapid>`_,
* IGS Ionosphere Working Group - `igsiono <ftp://igs-final.man.olsztyn.pl>`_.

The CDDIS has a large archive of IONEX data back to the 90's, but now requires authentication to access their data.
To use this service you must register a NASA EARTHDATA Account (https://urs.earthdata.nasa.gov/).
Then you need to create a ``~/.netrc`` file with the following text


`
machine urs.earthdata.nasa.gov login <username> password <password>
`


where ``<username>`` and ``<password>`` should be replaced with the appropriate values matching your account.
If using CDDIS, Spinifex will search for this file and raise an error if it does not exist.

On the servers, global ionospheric model data from different providers are stored, they can be identified with the
three letter prefix. The currently supported providers can be found in the how to download IONEX data example in:

 :doc:`examples/how_to_download`


We found that
for European purposes the high time resolution data from ``uqr`` give the best results, but this can be different for different
continents.
If your favorite server or data provider is not *yet* supported, please file an issue and we will try to include it.

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
* correct_uqrg_rms : [bool] The rms values for UQR maps are overestimated (see `Zhao et al. <https://link.springer.com/article/10.1007/s00190-021-01487-8>`_), use this flag to correct these. default=True


.. _tomion:

tomion
---------------------
The 2 layer ionospheric model provided by UPC-IonSat. This tomographic model is slower, but tends to give more stable results in long term pulsar observations.
See e.g. `Porayko, N. et al <https://link.springer.com/article/10.1007/s00190-023-01806-1>`_. Detailed information on the tomion fitting can be found
`here <_static/tomion_doc.pdf>`_. Note that from March 2022 the tomion files have high (15 minute) time resolution. Before that date only the 2 hour resolution data is
available, resulting in lower accuracy.

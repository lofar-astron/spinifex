# Spinifex

<!-- Re-enable when building these is working -->
<!-- ![Build status](git.astron.nl/spinifex/badges/main/pipeline.svg) -->
<!-- ![Test coverage](git.astron.nl/spinifex/badges/main/coverage.svg) -->

<!-- ![Latest release](https://git.astron.nl/templates/python-package/badges/main/release.svg) -->

<img src="_static/spinifex-logo.png" width="300" height="300" align="right" />

Pure Python tooling for ionospheric analyses in radio astronomy, e.g. getting total electron content and rotation
measure.

Spinifex is, in part, a re-write of [RMextract](https://github.com/lofar-astron/RMextract)
([Mevius, M. 2018](https://www.ascl.net/1806.024)). We have re-implemented all existing features of RMextract, but
`spinifex` is not directly backwards compatible with `RMextract`. We plan to extend to new features very soon.

Spinifex uses following external packages for RM modeling:

-   [PyIRI](https://doi.org/10.5281/zenodo.10139334) - Python implementation of the International Reference Ionosphere
    -   [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10139334.svg)](https://doi.org/10.5281/zenodo.10139334)
-   [ppigrf](https://github.com/IAGA-VMOD/ppigrf) - Pure Python IGRF (International Geomagnetic Reference Field)
    -   [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14231854.svg)](https://doi.org/10.5281/zenodo.14231854)

## Why 'Spinifex'?

[Spinifex](<https://en.wikipedia.org/wiki/Triodia_(plant)>) is a spiky grass native to the Australian continent. The
spines of spinifex are reminiscent of the ionospheric pierce points computed by this software. The 'spin' in spinifex
can also be thought to relate to the ionospheric Faraday rotation this software predicts.

## Command line interface

To run spinifex from the command line, use:

```
spinifex  --help
```

## Documentation

<div style="text-align:center"><img src="_static/altaz_example.png" width="300" height="300" />
<br />
<i>Example of RM output as function of azimuth and zenith angle</i>
</div>

Full documentation and examples of the Python module and the command-line tools are available on
[Read the Docs](https://spinifex.readthedocs.io/).

## Installation

Spinifex requires Python `>=3.9`. All dependencies can be installed with `pip`.

Stable release via PyPI:

```bash
pip install spinifex
```

If you have `casacore` installed, install with optional `python-casacore` dependency to work on MeasurementSets:

```bash
pip install spinifex[casacore]
```

Latest version on Gitlab:

```bash
pip install git+https://git.astron.nl/RD/spinifex
```

## Citing

If you use spinifex for your publications, please cite as:

> Mevius, M., Thomson, A., Dijkema, T.J.: Spinifex (2025),
> [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15000430.svg)](https://doi.org/10.5281/zenodo.15000430)

## License

This project is licensed under the Apache License Version 2.0

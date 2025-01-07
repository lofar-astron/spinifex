# Spinifex

<!-- Re-enable when building these is working -->
<!-- ![Build status](git.astron.nl/spinifex/badges/main/pipeline.svg) -->
<!-- ![Test coverage](git.astron.nl/spinifex/badges/main/coverage.svg) -->

<!-- ![Latest release](https://git.astron.nl/templates/python-package/badges/main/release.svg) -->

Pure Python tooling for ionospheric analyses in radio astronomy, e.g. getting total electron content and rotation
measure.

Spinifex is, in part, a re-write of [RMextract](https://github.com/lofar-astron/RMextract)
([Mevius, M. 2018](https://www.ascl.net/1806.024)). We have re-implemented all existing features of RMextract, but
`spinifex` is not directly backwards compatible with `RMextract`. We plan to extend to new features very soon.

Spinifex uses following reference ionosphere models:

-   [PyIRI](https://doi.org/10.5281/zenodo.10139334) - Python implementation of the International Reference Ionosphere
    -   [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10139334.svg)](https://doi.org/10.5281/zenodo.10139334)
-   [ppigrf](https://github.com/IAGA-VMOD/ppigrf) - Pure Python IGRF (International Geomagnetic Reference Field)
    -   [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14231854.svg)](https://doi.org/10.5281/zenodo.14231854)

### Why 'Spinifex'?

[Spinifex](<https://en.wikipedia.org/wiki/Triodia_(plant)>) is a spiky grass native to the Australian continent. The
spines of spinifex are reminiscent of the ionospheric pierce points computed by this software. The 'spin' in spinifex
can also be thought to relate to the ionospheric Faraday rotation this software predicts.

## Installation

Spinifex requires Python `>=3.9`. All dependencies can be installed with `pip`.

Stable release via PyPI:

```
pip install spinifex
```

Latest version on Gitlab:

```
pip install git+https://git.astron.nl/RD/spinifex
```

## Documentation

Full documentation and examples of the Python module and the command-line tools are available on
[Read the Docs](https://spinifex.readthedocs.io/).

## Contributing and Development

Test locally: `pip install tox`

With tox the same jobs as run on the CI/CD pipeline can be ran. These include unit tests and linting.

`tox`

To automatically apply most suggested linting changes execute:

`tox -e format`

## License

This project is licensed under the Apache License Version 2.0

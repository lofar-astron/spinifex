# Spinifex

![Build status](git.astron.nl/spinifex/badges/main/pipeline.svg)
![Test coverage](git.astron.nl/spinifex/badges/main/coverage.svg)

<!-- ![Latest release](https://git.astron.nl/templates/python-package/badges/main/release.svg) -->

Tool for ionospheric analyses, e.g. getting total electron content and rotation measure.

## Installation

```
pip install .
```

## Usage

```python
import spinifex
```

There are also command-line programs.

## Contributing

Test locally:
`pip install tox`

With tox the same jobs as run on the CI/CD pipeline can be ran. These include unit tests and linting.

`tox`

To automatically apply most suggested linting changes execute:

`tox -e format`

## License

This project is licensed under the Apache License Version 2.0

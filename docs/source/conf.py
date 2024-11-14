#  Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
#  SPDX-License-Identifier: Apache-2.0

import os

from spinifex import __version__

# -- General configuration ----------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.viewcode",
    "sphinxcontrib.apidoc",
    "sphinx_rtd_theme",
    "myst_parser"
]

# Assumes tox is used to call sphinx-build
project_root_directory = os.getcwd()

apidoc_module_dir = "../../spinifex"
apidoc_output_dir = "source_documentation"
apidoc_excluded_paths = []
apidoc_separate_modules = True
apidoc_toc_file = False
# This should include private methods but does not work
# https://github.com/sphinx-contrib/apidoc/issues/14
apidoc_extra_args = ["--private"]

# The suffix of source filenames.
source_suffix = [".rst"]

# The master toctree document.
master_doc = "index"

# General information about the project.
project = "spinifex"
copyright = "2023, ASTRON"

# openstackdocstheme options
repository_name = "git.astron.nl/spinifex"
bug_project = "none"
bug_tag = ""
html_last_updated_fmt = "%Y-%m-%d %H:%M"

# If true, '()' will be appended to :func: etc. cross-reference text.
add_function_parentheses = True

version = __version__

modindex_common_prefix = ["spinifex."]

# If true, the current module name will be prepended to all description
# unit titles (such as .. function::).
add_module_names = True

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = "sphinx"

# -- Options for HTML output --------------------------------------------------

# The theme to use for HTML and HTML Help pages.  Major themes that come with
# Sphinx are currently 'default' and 'sphinxdoc'.
# html_theme_path = ["."]
html_theme = "sphinx_rtd_theme"
html_static_path = ["static"]
html_css_files = [
    "css/custom.css",
]

# Output file base name for HTML help builder.
htmlhelp_basename = "%sdoc" % project

# Conf.py variables exported to sphinx rst files access using |NAME|
variables_to_export = [
    "project",
    "copyright",
    "version",
]

# Write to rst_epilog to export `variables_to_export` extract using `locals()`
# https://www.sphinx-doc.org/en/master/usage/configuration.html#confval-rst_epilog
frozen_locals = dict(locals())
rst_epilog = "\n".join(
    map(
        lambda x: f".. |{x}| replace:: {frozen_locals[x]}",  # noqa: F821
        variables_to_export,
    )
)
# Pep is not able to determine that frozen_locals always exists so noqa
del frozen_locals

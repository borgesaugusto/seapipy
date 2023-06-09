# Configuration file for the Sphinx documentation builder.
import sys
import os
project = 'SeapiPy'
author = 'Augusto Borges'

# release = '0.1'
version = '0.1.0'

# -- General configuration

sys.path.insert(0, os.path.abspath('../'))

extensions = [
    # 'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
]

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

# templates_path = ['_templates']

# -- Options for HTML output

html_theme = 'sphinx_rtd_theme'

# -- Options for EPUB output
# epub_show_urls = 'footnote'

autodoc_default_flags = ['members']
autosummary_generate = True

autodoc_mock_imports = ['numpy', 'scipy']
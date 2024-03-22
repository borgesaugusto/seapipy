# Configuration file for the Sphinx documentation builder.
import sys
import os
project = 'SeapiPy'
author = 'Augusto Borges'

# release = '0.1'
version = '0.2.1'

# -- General configuration

templates_path = ['_templates']
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


# -- Options for HTML output

html_theme = 'sphinx_rtd_theme'

# -- Options for EPUB output
# epub_show_urls = 'footnote'
autodoc_default_flags = ['members']
autosummary_generate = True
autosummary_imported_members = False
autoclass_content = "both"

autodoc_mock_imports = ['numpy', 'scipy', 'matplotlib']

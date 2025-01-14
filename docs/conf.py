# Configuration file for the Sphinx documentation builder.

# -- Project information

project = 'FaST'
copyright = '2024, Valerio Fulci'
author = 'Valerio Fulci'

release = '0.3'
version = '0.3.0'

# -- General configuration

extensions = [
    'sphinx.ext.duration',
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

templates_path = ['_templates']

# -- Options for HTML output


# -- Options for EPUB output
epub_show_urls = 'footnote'

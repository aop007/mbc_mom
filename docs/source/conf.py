# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'mbc-mom'
copyright = '2026, Alexis Ouellet Patenaude, ing.'
author = 'Alexis Ouellet Patenaude, ing.'
release = '0.1.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'matplotlib.sphinxext.plot_directive'
]

templates_path = ['_templates']
exclude_patterns = []

# Ensure autodoc looks at the type hints in your .pyi files
autodoc_typehints = "description"

plot_include_source = True
plot_html_show_source_link = False
plot_html_show_formats = False

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'pydata_sphinx_theme'
html_static_path = ['_static']

# Minimal Sphinx configuration for ReadTheDocs
# This is a dummy file required by ReadTheDocs
# Actual documentation is built using Julia's Documenter.jl

project = 'HyQMOM.jl'
copyright = '2024, Computational Physics Group'
author = 'Computational Physics Group'

# The master toctree document
master_doc = 'index'

# ReadTheDocs theme (not actually used since we build with Documenter.jl)
html_theme = 'alabaster'

# Disable most Sphinx extensions since we're not using Sphinx
extensions = []


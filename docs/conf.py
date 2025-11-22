# Minimal Sphinx configuration for ReadTheDocs
# This is required by ReadTheDocs but not used
# Actual documentation is built using Julia's Documenter.jl in post_build

project = 'HyQMOM.jl'
copyright = '2024, Computational Physics Group'
author = 'Computational Physics Group'

# Point to a non-existent master doc so Sphinx has nothing to build
master_doc = 'nonexistent'

html_theme = 'alabaster'
extensions = []

# Exclude all patterns to prevent Sphinx from processing anything
exclude_patterns = ['*']


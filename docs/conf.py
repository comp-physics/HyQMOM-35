# Minimal Sphinx configuration for ReadTheDocs
# This is required by ReadTheDocs but not actually used
# Actual documentation is built using Julia's Documenter.jl in post_build
# The post_build job will overwrite this Sphinx output with Documenter.jl HTML

project = 'HyQMOM.jl'
copyright = '2024, Computational Physics Group'
author = 'Computational Physics Group'

# Minimal master doc (will be overwritten by Documenter.jl)
master_doc = 'nonexistent'

html_theme = 'alabaster'
extensions = []

# Exclude markdown and Julia files from Sphinx
exclude_patterns = ['*.md', 'src/**', 'make.jl', 'Project.toml', 'Manifest.toml', 'build/**']


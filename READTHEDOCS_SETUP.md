# ReadTheDocs Setup Guide for HyQMOM.jl

This repository is configured for automatic documentation deployment to ReadTheDocs using Julia's Documenter.jl.

## Current Configuration

### Files Involved

1. **`.readthedocs.yaml`** - Main ReadTheDocs configuration
   - Installs Julia 1.10.5
   - Sets up the documentation environment
   - Builds docs using `docs/make.jl`
   - Outputs to ReadTheDocs' expected location

2. **`docs/conf.py`** - Minimal Sphinx configuration (required by ReadTheDocs)
   - This is a dummy file; actual docs are built with Documenter.jl
   - Required for ReadTheDocs to recognize the project

3. **`docs/make.jl`** - Julia Documenter.jl build script
   - Automatically detects ReadTheDocs environment
   - Outputs to `$READTHEDOCS_OUTPUT/html` when on RTD
   - Skips GitHub Pages deployment when on RTD

4. **`.github/workflows/docs.yml`** - GitHub Actions workflow
   - Builds docs on every push to `main` or `docs` branches
   - Creates artifacts for preview

## Setting Up ReadTheDocs (First Time)

### 1. Connect Repository to ReadTheDocs

1. Go to [https://readthedocs.org](https://readthedocs.org)
2. Sign in with your GitHub account
3. Click "Import a Project"
4. Select `comp-physics/HyQMOM.jl` from your repositories
5. Click "Next"

### 2. Configure Project Settings

In the ReadTheDocs project dashboard:

1. **Name**: HyQMOM.jl
2. **Repository URL**: `https://github.com/comp-physics/HyQMOM.jl`
3. **Default branch**: `main`
4. **Default version**: `latest`

### 3. Enable Builds

1. Go to "Admin" → "Advanced Settings"
2. Ensure these settings:
   - **Default branch**: `main`
   - **Privacy level**: Public (or as needed)
   - **Build pull requests**: ✓ (optional, useful for previews)

### 4. Trigger First Build

1. Go to "Builds" tab
2. Click "Build Version: latest"
3. Monitor the build log for any errors

## How It Works

### Build Process

When you push to `main` or trigger a build on ReadTheDocs:

1. **RTD clones the repo** and checks out the specified branch
2. **Installs Julia** (1.10.5) via tarball download
3. **Sets up docs environment**:
   ```bash
   julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=".")); Pkg.instantiate()'
   ```
4. **Builds documentation**:
   ```bash
   julia --project=docs docs/make.jl
   ```
5. **Documenter.jl** detects `READTHEDOCS_OUTPUT` and builds HTML to the correct location
6. **RTD serves** the generated HTML at your project's URL

### Environment Variables

The build automatically sets:
- `HYQMOM_SKIP_PLOTTING=true` - Skips heavy plotting dependencies
- `JULIA_NUM_THREADS=2` - Limits threads for faster builds
- `READTHEDOCS_OUTPUT` - RTD provides this; Documenter.jl uses it

### Dual Deployment

This repo supports both:
- **ReadTheDocs**: Automatic builds from main branch
- **GitHub Pages**: Via manual deployment or GitHub Actions

The `docs/make.jl` file automatically:
- Deploys to GitHub Pages when NOT on ReadTheDocs
- Skips GitHub Pages deployment when on ReadTheDocs

## Updating Documentation

### Local Build

Build docs locally to test changes:

```bash
# From repository root
./build_docs.sh

# Or manually
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=".")); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

View the built docs:
```bash
./serve_docs.sh
# Opens at http://localhost:8000
```

### Push to Deploy

1. Make changes to documentation in `docs/src/`
2. Commit and push to `main` branch
3. ReadTheDocs automatically rebuilds (usually within 2-5 minutes)
4. Check build status at https://readthedocs.org/projects/hyqmom-jl/builds/

## Build Logs and Debugging

### Viewing Build Logs

1. Go to your ReadTheDocs project dashboard
2. Click "Builds" tab
3. Click on a specific build to view detailed logs

### Common Issues

#### Julia Download Fails
- **Symptom**: Build fails during Julia installation
- **Fix**: Update Julia URL in `.readthedocs.yaml` to latest stable version

#### Package Dependencies Fail
- **Symptom**: `Pkg.instantiate()` errors
- **Fix**: Ensure `docs/Project.toml` and `docs/Manifest.toml` are up to date

#### Documenter.jl Build Fails
- **Symptom**: Errors during `docs/make.jl`
- **Fix**: Check for missing docstrings, broken links, or invalid syntax in `.md` files

#### Wrong Output Location
- **Symptom**: RTD can't find the HTML files
- **Fix**: Verify `docs/make.jl` uses `get(ENV, "READTHEDOCS_OUTPUT", ...)` correctly

## Webhooks and Automation

ReadTheDocs should automatically set up a webhook on your GitHub repository. To verify:

1. Go to GitHub repo → Settings → Webhooks
2. You should see a webhook pointing to `readthedocs.org`
3. It should trigger on push events

If missing, ReadTheDocs can regenerate it from the project's "Admin" → "Integrations" page.

## Documentation URLs

Once set up, your docs will be available at:

- **Latest (main branch)**: `https://hyqmom-jl.readthedocs.io/en/latest/`
- **Stable (tagged versions)**: `https://hyqmom-jl.readthedocs.io/en/stable/`
- **Specific versions**: `https://hyqmom-jl.readthedocs.io/en/v1.0.0/` (if tagged)

## Version Control

To build docs for specific versions:

1. Create a git tag: `git tag v1.0.0`
2. Push the tag: `git push origin v1.0.0`
3. ReadTheDocs will automatically build that version
4. Set it as "stable" in RTD admin panel if desired

## Troubleshooting

### Build Takes Too Long
- Julia download and package precompilation can take 5-10 minutes
- This is normal for Julia projects on ReadTheDocs

### PlotlyJS or GLMakie Errors
- Should be avoided by `HYQMOM_SKIP_PLOTTING=true`
- If not, check `docs/Project.toml` for plotting dependencies

### Links to GitHub Source
- Should work automatically via `repo` setting in `docs/make.jl`
- URL: `https://github.com/comp-physics/HyQMOM.jl/blob/{commit}{path}#{line}`

## Support

- **ReadTheDocs Docs**: https://docs.readthedocs.io/
- **Documenter.jl Docs**: https://documenter.juliadocs.org/
- **Issues**: Open an issue on the GitHub repository

---

**Last Updated**: November 2024 (after repository restructure)


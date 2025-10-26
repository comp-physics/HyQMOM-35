# Ensure heavy plotting deps are skipped during doc build
ENV["HYQMOM_SKIP_PLOTTING"] = get(ENV, "HYQMOM_SKIP_PLOTTING", "true")

using Pkg
# Make sure the docs env sees the local package
Pkg.develop(PackageSpec(path=joinpath(@__DIR__, "..")))
Pkg.instantiate()

using Documenter
using DocStringExtensions
using Literate
using HyQMOM

# Convert selected example(s) to tutorials without executing them (fast, CI-friendly)
# Temporarily disabled due to parsing issues with the example file
# tutorials_dir = joinpath(@__DIR__, "src", "tutorials")
# mkpath(tutorials_dir)
# examples_root = joinpath(@__DIR__, "..", "examples")
# for ex in ["run_3d_jets_timeseries.jl"]
#     src = joinpath(examples_root, ex)
#     if isfile(src)
#         Literate.markdown(src, tutorials_dir; execute=false, flavor=Literate.DocumenterFlavor())
#     end
# end

# Doctest setup
DocMeta.setdocmeta!(HyQMOM, :DocTestSetup, :(using HyQMOM); recursive=true)

# Determine the HTML output dir (Read the Docs provides READTHEDOCS_OUTPUT)
html_out = joinpath(get(ENV, "READTHEDOCS_OUTPUT", joinpath(@__DIR__, "build")), "html")

makedocs(
    sitename = "HyQMOM.jl",
    modules = [HyQMOM],
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", "false") == "true"),
    clean = true,
    pages = [
        "Home" => "index.md",
        "Quickstart" => "quickstart.md",
        "User Guide" => "user_guide.md",
        "MPI & Parallelization" => "mpi.md",
        "Tutorials" => [
            "Interactive Visualization" => "tutorials/interactive_visualization.md",
        ],
        "API Reference" => "api.md",
        "Developer Guide" => "dev_guide.md",
    ],
    build = html_out,
)

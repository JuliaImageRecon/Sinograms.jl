# push!(LOAD_PATH,"../src/")

using Sinograms
using Documenter

DocMeta.setdocmeta!(Sinograms, :DocTestSetup, :(using Sinograms); recursive=true)

makedocs(;
    modules = [Sinograms],
    authors = "Jeff Fessler <fessler@umich.edu> and contributors",
    repo = "https://github.com/JeffFessler/Sinograms.jl/blob/{commit}{path}#{line}",
    sitename = "Sinograms.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
#       canonical = "https://JeffFessler.github.io/Sinograms.jl/stable",
#       assets = String[],
    ),
    pages = [
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo = "github.com/JeffFessler/Sinograms.jl.git",
    devbranch = "main",
    devurl = "dev",
    versions = ["stable" => "v^", "dev" => "dev"]
#   push_preview = true,
)

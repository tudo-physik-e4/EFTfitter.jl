# Use
#
#     DOCUMENTER_DEBUG=true julia --color=yes make.jl local [nonstrict] [fixdoctests]
#
# for local builds.

using Documenter
using EFTfitter

makedocs(
    sitename = "EFTfitter.jl",
    modules = [EFTfitter],
    format = Documenter.HTML(
        prettyurls = !("local" in ARGS),
        canonical = "https://github.com/Cornelius-G/EFTfitter.jl"
    ),
    pages = [
        "Home" => "index.md",
        "Installation" => "installation.md",
        "Tutorial" => "tutorial.md",
        "Advanced Tutorial" => "advanced_tutorial.md",
        "BLUE Example" => "BLUE.md",
        "Plotting" => "plotting.md",
        "API" => "api.md",
        "LICENSE" => "LICENSE.md",
    ],
    doctest = ("fixdoctests" in ARGS) ? :fix : true,
    linkcheck = ("linkcheck" in ARGS),
    strict = !("nonstrict" in ARGS),
)

deploydocs(
    repo = "github.com/Cornelius-G/EFTfitter.jl.git",
    forcepush = true,
    push_preview = true,
    devbranch = "dev"
)

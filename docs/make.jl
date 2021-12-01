using DelaySSAdocs
using Documenter

DocMeta.setdocmeta!(DelaySSAdocs, :DocTestSetup, :(using DelaySSAdocs); recursive=true)

makedocs(;
    modules=[DelaySSAdocs],
    authors="Xiaoming Fu",
    repo="https://github.com/palmtree2013/DelaySSAdocs.jl/blob/{commit}{path}#{line}",
    sitename="DelaySSAdocs.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://palmtree2013.github.io/DelaySSAdocs.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/palmtree2013/DelaySSAdocs.jl",
    devbranch="main",
)

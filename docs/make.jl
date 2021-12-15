using DelaySSAdocs
using Documenter

DocMeta.setdocmeta!(DelaySSAdocs, :DocTestSetup, :(using DelaySSAdocs); recursive=true)

makedocs(;
    modules=[DelaySSAdocs],
    authors="Xiaoming Fu",
    repo="https://github.com/palmtree2013/DelaySSAdocs.jl/blob/{commit}{path}#{line}",
    sitename="DelaySSAdocs.jl",
    format=Documenter.HTML(;
        mathengine=Documenter.Writers.HTMLWriter.MathJax2(),
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://palmtree2013.github.io/DelaySSAdocs.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Tutorials" => [
            "tutorials/tutorials.md", 
            "tutorials/bursty.md",
            "tutorials/delay_degradation.md",
            # "tutorials/delay_multidegradation.md",
            "tutorials/heterogeneous_delay.md",
            "tutorials/delay_oscillator.md",
            "tutorials/stochastic_delay.md",
        ],
        "Algorithm" => [
            "algorithms/notations.md",      
            "algorithms/delayrejection.md",      
            "algorithms/delaydirect.md",
            "algorithms/delaymnrm.md", 
        ],
        "Theory" => "theory.md",
        "API" => "api.md",
    ],
)

deploydocs(;
    repo="github.com/palmtree2013/DelaySSAdocs.jl",
    devbranch="main",
)

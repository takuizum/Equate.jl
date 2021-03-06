using Documenter
using Equate

makedocs(
    sitename = "Equate.jl",
    format = Documenter.HTML(
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://takuizum.github.io/Equate.jl/stable",
        assets=String[],
    ),
    modules = [Equate],
    pages = ["Home" => "index.md"]
)


deploydocs(;
    repo="github.com/takuizum/Equate.jl",
)
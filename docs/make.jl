using ControlSystemsTools
using Documenter

DocMeta.setdocmeta!(ControlSystemsTools, :DocTestSetup, :(using ControlSystemsTools); recursive=true)

makedocs(;
    modules=[ControlSystemsTools],
    authors="F. Capolupo",
    sitename="ControlSystemsTools.jl",
    format=Documenter.HTML(;
        canonical="https://FraCpl.github.io/ControlSystemsTools.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/FraCpl/ControlSystemsTools.jl",
    devbranch="master",
)

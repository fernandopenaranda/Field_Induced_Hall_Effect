using Field_Induced_Hall_Effect
using Documenter

DocMeta.setdocmeta!(Field_Induced_Hall_Effect, :DocTestSetup, :(using Field_Induced_Hall_Effect); recursive=true)

makedocs(;
    modules=[Field_Induced_Hall_Effect],
    authors="Fernando Penaranda <fernandopenaranda@github.com> and contributors",
    sitename="Field_Induced_Hall_Effect.jl",
    format=Documenter.HTML(;
        canonical="https://fernandopenaranda.github.io/Field_Induced_Hall_Effect.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/fernandopenaranda/Field_Induced_Hall_Effect.jl",
    devbranch="main",
)

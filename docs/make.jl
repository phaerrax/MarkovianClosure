using MarkovianClosure
using Documenter

DocMeta.setdocmeta!(MarkovianClosure, :DocTestSetup, :(using MarkovianClosure); recursive=true)

makedocs(;
    modules=[MarkovianClosure],
    authors="Davide Ferracin <davide.ferracin@protonmail.com> and contributors",
    sitename="MarkovianClosure.jl",
    format=Documenter.HTML(;
        canonical="https://phaerrax.github.io/MarkovianClosure.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/phaerrax/MarkovianClosure.jl",
    devbranch="main",
)

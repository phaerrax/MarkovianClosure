using Documenter
using ITensors
using MarkovianClosure

DocMeta.setdocmeta!(
    MarkovianClosure, :DocTestSetup, :(using MarkovianClosure); recursive=true
)

makedocs(;
    modules=[MarkovianClosure],
    checkdocs=:exported,
    authors="Davide Ferracin <davide.ferracin@protonmail.com> and contributors",
    sitename="MarkovianClosure.jl",
    format=Documenter.HTML(;
        canonical="https://phaerrax.github.io/MarkovianClosure.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=["Home" => "index.md"],
)

deploydocs(;
    branch="gh-pages", repo="github.com/phaerrax/MarkovianClosure.jl", devbranch="main"
)

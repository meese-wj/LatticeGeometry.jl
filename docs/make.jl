using LatticeGeometry
using Documenter

DocMeta.setdocmeta!(LatticeGeometry, :DocTestSetup, :(using LatticeGeometry); recursive=true)

makedocs(;
    modules=[LatticeGeometry],
    authors="W. Joe Meese <meesewj@gmail.com> and contributors",
    sitename="LatticeGeometry.jl",
    format=Documenter.HTML(;
        canonical="https://meese-wj.github.io/LatticeGeometry.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/meese-wj/LatticeGeometry.jl",
    devbranch="main",
)

using AllenNeuropixelsBase
using Documenter

DocMeta.setdocmeta!(AllenNeuropixelsBase, :DocTestSetup, :(using AllenNeuropixelsBase); recursive=true)

makedocs(;
    modules=[AllenNeuropixelsBase],
    authors="brendanjohnharris <brendanjohnharris@gmail.com> and contributors",
    repo="https://github.com/brendanjohnharris/AllenNeuropixelsBase.jl/blob/{commit}{path}#{line}",
    sitename="AllenNeuropixelsBase.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(
    repo = "github.com/brendanjohnharris/AllenNeuropixelsBase.jl.git",
)

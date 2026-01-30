using CLOUD18_calibration
using Documenter

DocMeta.setdocmeta!(CLOUD18_calibration, :DocTestSetup, :(using CLOUD18_calibration); recursive=true)

makedocs(;
    modules=[CLOUD18_calibration],
    authors="Timo Wittler wittler.timo@uibk.ac.at, Clea Ruth clea.ruth@uibk.ac.at",
    sitename="CLOUD18_calibration.jl",
    format=Documenter.HTML(;
        canonical="https://wimolimo.github.io/CLOUD18_calibration.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/wimolimo/CLOUD18_calibration.jl",
    devbranch="master",
)

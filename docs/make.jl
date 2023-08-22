push!(LOAD_PATH,"../src/")
using SeismicQ
using Documenter
makedocs(
         sitename = "SeismicQ.jl",
         modules  = [SeismicQ],
         pages=[
                "Home" => "index.md"
               ])
deploydocs(;
    repo="github.com/tduretz/SeismicQ.jl",
)
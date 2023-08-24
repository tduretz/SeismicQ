module SeismicQ

include("Sources.jl")
export Ricker

include("Spec.jl")
export Getfreq,Spec,ComputeQgraph

include("GenFakeData.jl")
export GenMatrix

end # module SeismicQ

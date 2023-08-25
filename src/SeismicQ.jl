module SeismicQ

include("Sources.jl")
export Ricker

include("Rheology.jl")
export f_bulk
export f_shear
export f_relax
include("Spec.jl")
export Getfreq, Spec, ComputeQgraph

include("GenFakeData.jl")
export GenMatrix

end # module SeismicQ

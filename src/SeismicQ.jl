module SeismicQ

include("Sources.jl")
export Ricker

include("Rheology.jl")
export f_bulk
export f_visc
export f_shear
export f_relax
include("Spec.jl")
export Getfreq, Spec, ComputeQgraph

include("GenFakeDataCleanFinal.jl")
export GenAttenuatedRicker

include("PlotReceiverGather.jl")
export PlotReceiverGather

end # module SeismicQ

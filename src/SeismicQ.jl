module SeismicQ

include("Sources.jl")
export Ricker

include("Rheology.jl")
export θs, χs, θb, χb

include("Spec.jl")
export Getfreq, Spec, ComputeQgraph

include("GenFakeDataCleanFinal.jl")
export GenAttenuatedRicker

include("PlotReceiverGather.jl")
export PlotReceiverGather

end # module SeismicQ

module SeismicQ

include("Sources.jl")
export Ricker

include("Spec.jl")
export Getfreq,Spec,ComputeQgraph

include("GenFakeDataCleanFinal.jl")
export GenAttenuatedRicker

include("PlotReceiverGather.jl")
export PlotReceiverGather

end # module SeismicQ

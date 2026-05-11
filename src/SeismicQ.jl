module SeismicQ

include("Sources.jl")
export Ricker,ExplosiveX,ExplosiveY,Gaussian, 
       DoubleCoupleX, DoubleCoupleY, RotationalX, RotationalY   

include("Rheology.jl")
export θs, χs, ηs, θb, χb, ηb

include("Spec.jl")
export Getfreq, Spec, ComputeQgraph

include("GenFakeDataCleanFinal.jl")
export GenAttenuatedRicker

include("PlotReceiverGather.jl")
export PlotReceiverGather

include("io_sismo.jl")
export make_filename
export write_sismo_file
export read_sismo_file

end # module SeismicQ

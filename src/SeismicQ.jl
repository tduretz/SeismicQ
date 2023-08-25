module SeismicQ

include("Sources.jl")
export Ricker

include("Rheology.jl")
export f_bulk
export f_shear
export f_relax

end # module SeismicQ

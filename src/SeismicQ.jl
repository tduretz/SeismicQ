module SeismicQ

include("Sources.jl")
export Ricker

include("spec.jl")
export getfreq,spec,compute_qgraph

include("gen_fake_data.jl")
export gen_matrix

end # module SeismicQ

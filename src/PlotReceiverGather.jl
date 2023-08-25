using SeismicQ, Plots, SpecialFunctions, LinearAlgebra, Printf

@doc raw"""
    PlotReceiverGather(listₓ,time_vec,acc_vec) 

Function that plots a receiver gather 
The function takes as inputs:
        listₓ : vector of geophone distances to the source [m]
        time_vec : vector containing the time steps of the wave signal [s]
        acc_vec : matrix of wave acceleration at each geophone position

# Examples
```julia-repl
julia>  PlotReceiverGather(0:1000:5000,time_vec,acc_vec) 

```
"""
function PlotReceiverGather(listₓ,time_axis,acc_vec) 

    p = heatmap(listₓ,time_axis,acc_vec', 
    color=palette(:RdBu, 100, rev=true),
    clim=(-0.5, 0.5),
    cbar=true,
    label=" ", 
    yflip=true,
    title = "Receiver gather", 
    xlabel = "Geophone distance to the source [m]",
    ylabel = "Time [s]")
    display(plot(p))
end
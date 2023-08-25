using SeismicQ

function main()

    # Geophone position [m]
    listₓ = 0:100:5000;

    # Time domain
    Δt  = 1e-3
    Nt  = 2000

    # Central frequency of the source [Hz]
    𝑓₀  = 10.

    # Velocities
    Vp = 7000 # m/s
    Vs = 4000 
    αp = 2e-4
    αs = 4e-4 # (π * f)/ (Q * V)

    time_vec,acc_vec=GenAttenuatedRicker(listₓ,Δt,Nt,Vp,Vs,αp,αs,𝑓₀)

    PlotReceiverGather(listₓ,time_vec,acc_vec)
end

main()

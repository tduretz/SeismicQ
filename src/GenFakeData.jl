# Function that generates a seismic trace (Ricker function), attenuating along time, 
# and plot a receiver gather (received waves at different geophone positions)

using SeismicQ, Plots, SpecialFunctions, LinearAlgebra, Printf

function GenMatrix(Vp,Vs)

    # Functtion that create a vector of time and 
    function FausseTrace(x,Δt,Nt,t,Vp,Vs,αp,αs,𝑓₀,t₀)

        tₓp = x/Vp
        tₓs = x/Vs

        # Storage
        time = zeros(Nt)
        acc = zeros(Nt)

        # Time loop
        for it=1:Nt

            # Compute Ricker function
            t += Δt
            a = 0.5 * exp(-αp*x)*Ricker(t, tₓp+t₀, 𝑓₀)+ 0.5 * exp(-αs*x)*Ricker(t, tₓs+t₀, 𝑓₀)

            # For visualisation purpose
            time[it] = t
            acc[it]  = a
        end

        return(time',acc')

    end

    # Geophone position [m]
    listₓ = 0:100:5000;

    # Time domain
    Δt  = 1e-3
    Nt  = 2000
    # Central frequency of the source [Hz]
    𝑓₀  = 10.
    t₀ = 1.0/𝑓₀
    t   = -t₀

    # Velocities
    #Vp = 7000 # m/s
    #Vs = 4000 
    αp = 2e-4
    αs = 4e-4 # (π * f)/ (Q * V)

    time_axis = t:(Δt*Nt-t)/(Nt-1):Δt*Nt ;
    time_vec = zeros(size(listₓ,1),Nt);
    acc_vec  = zeros(size(listₓ,1),Nt);
    time_vec[1,:],acc_vec[1,:]= FausseTrace(listₓ[1,1],Δt,Nt,t,Vp,Vs,αp,αs,𝑓₀,t₀);
    dist_axis = 0:1: size(listₓ,1)-1; # distances des geophones pour l instant 1,2,3,...

    for i=1:size(listₓ,1)-1
        time_vec[i+1,:],acc_vec[i+1,:] = FausseTrace(listₓ[i+1,1],Δt,Nt,t,Vp,Vs,αp,αs,𝑓₀,t₀)
    end

    seismic_matrix = hcat(listₓ,acc_vec);

    return seismic_matrix
end
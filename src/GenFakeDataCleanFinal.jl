# Function that generates a matrix that corresponds to the seismic trace (as a Ricker function), attenuated along time
# and received at different geophone distances to the source
# Inputs :
    # listₓ : geophone distances to the source (vector) e.g. listₓ = 0:100:5000 ;
    # Δt : time step of the wave [s] e.g. Δt  = 1e-3 ;
    # Nt : number of time steps of the wave e.g. Nt  = 2000 ;
    # Vp : P-wave velocity [m/s] e.g. Vp = 7000 ;
    # Vs : S-wave velocity [m/s] e.g. Vp = 4000 ;
    # αp : attenuation factor for P-wave e.g. αp = 2e-4 ;
    # αs : attenuation factor for S-wave e.g. αs = 4e-4 ;
    # 𝑓₀ : central frequency of the source e.g. 𝑓₀  = 10.0 ; 
# Outputs : 
    # time_axis : time points of the received waves (vector)
    # acc_vec : matrix of wave acceleration at each time points and the geophone positions 

using SeismicQ, Plots, SpecialFunctions, LinearAlgebra, Printf

@doc raw"""
    time_vec,acc_vec = GenAttenuatedRicker(listₓ,Δt,Nt,Vp,Vs,αp,αs,𝑓₀)  ;

Function that generates a seismic trace (Ricker wavelet) of P and S waves.
This wave is attenuated along time and the function creates a matrix of wave amplitude as a function of time and distance to the source. 
The function takes as inputs:
    listₓ : vector of geophone distances to the source [m]
    Δt : time step of the wave signal [s] 
    Nt : number of time steps of the wave signal 
    Vp : P-wave velocity [m/s] 
    Vs : S-wave velocity [m/s] 
    αp : attenuation factor for P-wave 
    αs : attenuation factor for S-wave 
    𝑓₀ : central frequency of the source [Hz]  

and return:
    time_vec : vector containing the time steps of the received wave signal
    acc_vec  : matrix of wave acceleration at each geophone position

# Examples
```julia-repl
julia>  time_vec,acc_vec = GenAttenuatedRicker(0:1000:5000,1e-3,2000,7000,4000,2e-4,4e-4,10.0) 

```
"""
function GenAttenuatedRicker(listₓ,Δt,Nt,Vp,Vs,αp,αs,𝑓₀) 

    # Central period of the source [Hz]
    t₀ = 1.0/𝑓₀ ;
    t   = -t₀ ;

    # Vectors initialization
    time_vec  = zeros(Nt);
    acc_vec   = zeros(size(listₓ,1),Nt);
    tₓp       = zeros(size(listₓ,1),1);
    tₓs       = zeros(size(listₓ,1),1);

    # Creation of Ricker traces and attenuated traces at geophone locations specified by listₓ
    for i=1:size(listₓ,1)

        # P and S wave times as a function of their imput velocities
        tₓp[i,1] = listₓ[i,1]/Vp ;
        tₓs[i,1] = listₓ[i,1]/Vs ;
        t = -t₀

            # Time loop
            for it=1:Nt

                # Compute Ricker function with attenuation factor for P (αp) and S (αs) waves
                t += Δt ;
                time_vec[it]   = t
                acc_vec[i,it]  = 0.5 * exp(-αp*listₓ[i,1])*Ricker(t, tₓp[i,1]+t₀, 𝑓₀)+ 
                          0.5 * exp(-αs*listₓ[i,1])*Ricker(t, tₓs[i,1]+t₀, 𝑓₀) ;
            end
    end

    # Outputs
    return time_vec,acc_vec

    
end


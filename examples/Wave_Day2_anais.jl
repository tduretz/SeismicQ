using SeismicQ, Plots

function MainSource()
    
    # Spatial extent
    Lx  = 25.0      # (en m)

    # Mechanical parameters

    ρ₀ = 1100.0      # (kg.m-3)
    K₀ = 1.0e9       # (Pa)
    G₀ = 1.0e8       # (Pa)



    # Discretization
    Ncx = 10
    Δx  = Lx/Ncx

    # Central frequency of the source [Hz]
    𝑓₀  = 10.
    t₀  = 1.2/𝑓₀

    # Time domain
    Δt  = 1e-3
    Nt  = 1000
    t   = -t₀
    v   = 0.0

    # Storage

    szv  = (Ncx+1,)
    szc  = (Ncx+2,)
    K    = ones(Ncx+2)*K₀
    G    = ones(Ncx+2)*G₀
    ρ    = ones(Ncx+1)*ρ₀
    Vx   = zeros(Ncx+1)
    ε̇    = (xx = zeros(szc), yy = zeros(szc) , zz = zeros(szc)  )
    τ    = ( xx = zeros(szc), yy = zeros(szc) , zz = zeros(szc) )
    P    = zeros(Ncx+2) 
    ∇V   = zeros(Ncx+2)
    
 #    time_src = zeros(Nt)
 #    acc_src  = zeros(Nt)
 #    vel_src  = zeros(Nt)
 #    v_src    = 0.


    @show ε̇.xx[1]

    # Time loop
    for it=1:Nt

        V0 = V 

        # Compute Ricker function
        t += Δt
        a  = Ricker(t, t₀, 𝑓₀)
        v += a*Δt
        
        # BC 
        V.x[1]  = v

        # gard V components

        ∂Vx∂x   = (V.x[2:end] - V.x[1:end-1]) / Δx
        ∂Vy∂y   = (V.y[2:end] - V.y[1:end-1]) / Δy
        ∂Vz∂z   = (V.z[2:end] - V.z[1:end-1]) / Δz

        ∂Vx∂y   = (V.y[2:end] - V.y[1:end-1]) / Δy
        ∂Vx∂z   = (V.y[2:end] - V.y[1:end-1]) / Δy

        ∂Vy∂z   = (V.y[2:end] - V.y[1:end-1]) / Δz
        ∂Vz∂y   = (V.z[2:end] - V.z[1:end-1]) / Δy


        # divergence 
        ∇V      = ∂Vx∂x +  ∂Vx∂y + ∂Vx∂z

        # Deviatoric strain rate 

         ε̇.xx[2:end-1] = ∂Vx∂x - 1/3 * ∇V[2:end-1]
         ε̇.yy[2:end-1] = ∂Vy∂y - 1/3 * ∇V[2:end-1]
         ε̇.zz[2:end-1] = ∂Vz∂z - 1/3 * ∇V[2:end-1] 
         ε̇.xy[2:end-1] = 0.5 * (∂Vx∂y + ∂Vy∂x)
         ε̇.xz[2:end-1] = 0.5 * (∂Vx∂z + ∂Vz∂x)
         ε̇.yz[2:end-1] = 0.5 * (∂Vy∂z + ∂Vz∂y)
   
        # Stress update 
       τ.xx = f_shear(G) * Δt * (ε̇.xx) + f_relax(G) * τ.xx
       τ.yy = f_shear(G) * Δt * (ε̇.yy) + f_relax(G) * τ.yy
       τ.zz = f_shear(G) * Δt * (ε̇.zz) + f_relax(G) * τ.zz
       τ.xy = f_shear(G) * Δt * (ε̇.xy) + f_relax(G) * τ.xy
       τ.xz = f_shear(G) * Δt * (ε̇.xz) + f_relax(G) * τ.xz
       τ.yz = f_shear(G) * Δt * (ε̇.yz) + f_relax(G) * τ.yz


        # Pressure update
        P   = P0 + Δt * f_bulk(K) * ∇V
        dvdt * rho = div τ - grad P 
        V.x = V0.x + Δt / ρ * ( (τ.xx[2:end]-τ.xx[1:end-1]) / Δx
                               +(τ.xy[2:end]-τ.xy[1:end-1]) / Δy
                               +(τ.xy[2:end]-τ.x[1:end-1]) / Δz
                               - P[2;end]/Δx
        
        V.y = V0.y + Δt / ρ * ( (τ.yy[2:end]-τ.yy[1:end-1]) / Δx
                               +(τ.xy[2:end]-τ.xy[1:end-1]) / Δy
                               +(τ.xy[2:end]-τ.x[1:end-1]) / Δz
                               - P[2;end]/Δx
        V.z = V0.z + Δt / ρ * ( (τ.xx[2:end]-τ.xx[1:end-1]) / Δx
                               +(τ.xy[2:end]-τ.xy[1:end-1]) / Δy
                               +(τ.xy[2:end]-τ.x[1:end-1]) / Δz
                               - P[2;end]/Δx



        # For visualisation purpose
 #       time[it] = t
 #       acc[it]  = a
 #       vel[it]  = v
    end

    # Visualisation
 #   p1 = plot(time, acc, xlabel="t", ylabel="a")
 #   p2 = plot(time, vel, xlabel="t", ylabel="v")
 #   plot(p1, p2, layout=(2,1))

end

function f_bulk(K)
    return Ṗ = -K
end

function f_shear(G)
    return 2*G 
end

function f_relax(G)
    return 1.
end


MainSource()
using SeismicQ, Plots

function MainSource()
    
    # Spatial extent
    Lx  = 250.0

    # Mechanical parameters 
    ρ₀   = 1500.0
    K₀   = 1.e9
    G₀   = 1.e8
     
    # Discretization
    Ncx = 100
    Δx  = Lx/Ncx
    Δy  = 1.
    Δz  = 1.
    xv  = LinRange(0,Lx,Ncx+1)
    # Source parameters
    # Central frequency of the source [Hz]
    𝑓₀  = 10.
    t₀  = 1.2/𝑓₀

    # Time domain
    Δt   = min(1e10, 2.2*sqrt(ρ₀/K₀)) # Courant criteria from wavespeed
    Nt   = 500
    Nout = 50
    t    = -t₀
    v    = 0.0
   
    # Storage on centers # +2 for ghost nodes for BCs
    szv = (Ncx+1,)
    szc = (Ncx+2,)
    #
    K    = ones(Ncx+2)*K₀ 
    # G
    G    = ones(Ncx+2)*G₀
    # ε̇xx
    ε̇ = ( xx=zeros(szc), yy=zeros(szc), zz=zeros(szc), xy=zeros(szc), yz=zeros(szc), xz=zeros(szc) )  
    # ∇V
    ∇V   = zeros(Ncx+2)
    # P 
    P   = zeros(Ncx+2)
    # τ
    τ = ( xx=zeros(szc), yy=zeros(szc), zz=zeros(szc), xy=zeros(szc), yz=zeros(szc), xz=zeros(szc) )  
    ∂Vx∂x = zeros(szc)
    # Storage on vertices
    # Vx
    V  = ( x=zeros(szv), y=zeros(szv), z=zeros(szv))
    # ρ is on vx
    ρ = ones(szv)*ρ₀ 

    
    # Time loop
    for it=1:Nt
        # V0 = V  # !!!!!!!!!!!!!!! MEGA ACHTUNG!
        # Compute Ricker function
        t     += Δt
        a      = Ricker(t, t₀, 𝑓₀)
        v     += a*Δt
        V.x[1] = v

        # grad V components
        @. ∂Vx∂x[2:end-1] = (V.x[2:end]-V.x[1:end-1])/Δx
        # ∂Vy∂y = (V.y[2:end]-V.y[1:end-1])/Δy
        # ∂Vz∂z = (V.z[2:end]-V.z[1:end-1])/Δz

        # ∂Vx∂y = (V.x[2:end]-V.x[1:end-1])/Δy
        # ∂Vy∂x = (V.y[2:end]-V.y[1:end-1])/Δx
        
        # ∂Vx∂z = (V.x[2:end]-V.x[1:end-1])/Δz
        # ∂Vz∂x = (V.z[2:end]-V.z[1:end-1])/Δx

        # ∂Vy∂z = (V.y[2:end]-V.y[1:end-1])/Δz
        # ∂Vz∂y = (V.z[2:end]-V.z[1:end-1])/Δy
        
        # divergence
        @. ∇V = ∂Vx∂x #+∂Vy∂y+∂Vz∂z
        # deviatoric strain rate 
        @. ε̇.xx = ∂Vx∂x - 1/3*∇V
        #  ε̇.yy[2:end-1] = ∂Vy∂y - 1/3*∇V[2:end-1] 
        #  ε̇.zz[2:end-1] = ∂Vz∂z - 1/3*∇V[2:end-1] 
        #  ε̇.xy[2:end-1] = 0.5*(∂Vx∂y+∂Vy∂x)
        #  ε̇.xz[2:end-1] = 0.5*(∂Vx∂z+∂Vz∂x) 
        #  ε̇.yz[2:end-1] = 0.5*(∂Vz∂y+∂Vy∂z)

        # Stress update
        @. τ.xx = f_shear(G)*Δt*(ε̇.xx) + f_relax(G)*τ.xx
        # τ.yy = f_shear(G)*Δt*(ε̇.yy) + f_relax(G)*τ.yy
        # τ.zz = f_shear(G)*Δt*(ε̇.zz) + f_relax(G)*τ.zz
        # τ.xy = f_shear(G)*Δt*(ε̇.xy) + f_relax(G)*τ.xy
        # τ.xz = f_shear(G)*Δt*(ε̇.xz) + f_relax(G)*τ.xz
        # τ.yz = f_shear(G)*Δt*(ε̇.yz) + f_relax(G)*τ.yz

        # Pressure update 
        @. P -= Δt*f_bulk(K)*∇V

        # And now sum of the forces are equal to mass times acceleration
        # @. V.x = V0.x + Δt/ρ*((τ.xx[2:end]-τ.xx[1:end-1])/Δx 
        #                     +(τ.xy[2:end]-τ.xy[1:end-1])/Δy
        #                     +(τ.xz[2:end]-τ.xz[1:end-1])/Δz
        #                     - (P[2:end]-P[1:end-1])/Δx)
        @. V.x[2:end-1] += Δt/ρ[2:end-1]*((τ.xx[3:end-1]-τ.xx[2:end-2])/Δx - (P[3:end-1]-P[2:end-2])/Δx)

        # V.y = V0.y + Δt/rho*((τ.xy[2:end]-τ.xy[1:end-1])/Δx 
        #                     +(τ.yy[2:end]-τ.yy[1:end-1])/Δy
        #                     +(τ.yz[2:end]-τ.yz[1:end-1])/Δz
        #                     - (P[2:end]-P[1:end])/Δy)  
        # V.z = V0.z + Δt/rho*((τ.xz[2:end]-τ.xz[1:end-1])/Δx 
        #                     +(τ.yz[2:end]-τ.yz[1:end-1])/Δy
        #                     +(τ.zz[2:end]-τ.zz[1:end-1])/Δz
        #                     - (P[2:end]-P[1:end])/Δz)                    

        if mod(it, Nout)==0
            display(plot(xv, V.x))
            sleep(0.1)
        end
    end
end

function f_bulk(K) 
   return K
end

function f_shear(G)
    return 2*G
end
function f_relax(G)
    return 1.
end

MainSource()
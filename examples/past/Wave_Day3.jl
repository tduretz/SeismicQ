using SeismicQ, Plots

function MainSource()
    
    # Spatial extent
    Lx  = 25.0

    # Mechanical parameters 
    ρ₀   = 1500.0
    K₀   = 1.e9
    G₀   = 1.e8
    c₀   = sqrt((K₀+4/3*G₀)/ρ₀) 
     
    # Discretization
    Ncx = 1000
    Δx  = Lx/Ncx
    xv  = LinRange(0,Lx,Ncx+1)

    # Source parameters
    𝑓₀  = 100     # Central frequency of the source [Hz]
    t₀  = 1.2/𝑓₀

    # Time domain
    Δt   = min(1e10, Δx/c₀) # Courant criteria from wavespeed
    Nt   = 5000
    Nout = 100
    t    = -t₀
    v  = 0.0
   
    # Storage on centers # +2 for ghost nodes for BCs
    szv   = (Ncx+1,)
    szc   = (Ncx+2,)
    # Storage on centroids 
    K     = ones(szc)*K₀ 
    G     = ones(szc)*G₀
    ε̇     = ( xx=zeros(szc), yy=zeros(szc), zz=zeros(szc), xy=zeros(szc), yz=zeros(szc), xz=zeros(szc) )  
    ∇V    = zeros(Ncx+2)
    P     = zeros(Ncx+2)
    τ     = ( xx=zeros(szc), yy=zeros(szc), zz=zeros(szc), xy=zeros(szc), yz=zeros(szc), xz=zeros(szc) )  
    ∂Vx∂x = zeros(szc)
    # Storage on vertices
    V     = ( x=zeros(szv), y=zeros(szv), z=zeros(szv))
    ρ     = ones(szv)*ρ₀ 

    # Time loop
     @time for it=1:Nt

        # Compute Ricker function
        t     += Δt
        a      = Ricker(t, t₀, 𝑓₀)
        v     += a*Δt
        V.x[1] = v
       
        # Laetitia is not yet absorbing
        # vbc    = V.x[end]-c₀/Δx*Δt*(V.x[end]-V.x[end-1])
        # V.x[end] = vbc
        # @show vbc

        # Velocity gradient components
        @. ∂Vx∂x[2:end-1] = (V.x[2:end] - V.x[1:end-1])/Δx
        
        # Divergence
        @. ∇V   = ∂Vx∂x

        # Deviatoric strain rate 
        @. ε̇.xx = ∂Vx∂x - 1/3*∇V
      
        # Stress update
        @. τ.xx = f_shear(G)*Δt*(ε̇.xx) + f_relax(G)*τ.xx

        # Pressure update 
        @. P    = P - Δt*f_bulk(K)*∇V

        # Linear momentum balance
        @. V.x[2:end-1] = V.x[2:end-1] + Δt/ρ[2:end-1]*((τ.xx[3:end-1]-τ.xx[2:end-2])/Δx - (P[3:end-1]-P[2:end-2])/Δx)

        # Visualisation
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
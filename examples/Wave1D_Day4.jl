using SeismicQ, Plots

function MainSource()
    
    # Spatial extent
    Lx   = 100.0

    # Source parameters
    𝑓₀   = 200     # Central frequency of the source [Hz]
    t₀   = 1.2/𝑓₀
    isrc = 2 #Int((Ncx/2)+1)

    # Visualization
    visu = true

    # Absorbing boundaries
    cerW = false
    cerE = false

    # Mechanical parameters 
    ρ₀   = 1500.0
    K₀   = 1.e9
    G₀   = 1.e8
    c₀   = sqrt((K₀+4/3*G₀)/ρ₀) 
    De   = 1.0 # Deborah number
    η₀   = De*G₀ / 𝑓₀
     
    # Discretization
    Ncx = 500
    Δx  = Lx/Ncx
    xv  = LinRange(0,Lx,Ncx+1)
    xc  = LinRange(0-Δx/2,Lx+Δx/2,Ncx+2)

    # Time domain
    Δt   = min(1e10, Δx/c₀/2.1/1) # Courant criteria from wavespeed
    Nt   = 2000
    Nout = 50
    t    = -t₀
   
    # Storage on centers # +2 for ghost nodes for BCs
    szv   = (Ncx+1,)
    szc   = (Ncx+2,)
    # Storage on centroids 
    K     = ones(szc)*K₀ 
    G     = ones(szc)*G₀
    η     = ones(szc)*η₀
    ε̇     = ( xx=zeros(szc), yy=zeros(szc), zz=zeros(szc), xy=zeros(szc), yz=zeros(szc), xz=zeros(szc) )  
    ∇V    = zeros(Ncx+2)
    P     = zeros(Ncx+2)
    τ     = ( xx=zeros(szc), yy=zeros(szc), zz=zeros(szc), xy=zeros(szc), yz=zeros(szc), xz=zeros(szc) )  
    ∂Vx∂x = zeros(szc)
    # Storage on vertices
    V     = ( x=zeros(szv), y=zeros(szv), z=zeros(szv))
    ρ     = ones(szv)*ρ₀ 
    f_ext = zeros(szv)

    # BC
    Lbc        = 2.
    bc_filtW_v = 1.0 .- cerW.*exp.(-(xv.-0Lx).^2/Lbc.^2)
    bc_filtW_c = 1.0 .- cerW.*exp.(-(xc.-0Lx).^2/Lbc.^2)
    bc_filtE_v = 1.0 .- cerE.*exp.(-(xv.- Lx).^2/Lbc.^2)
    bc_filtE_c = 1.0 .- cerE.*exp.(-(xc.- Lx).^2/Lbc.^2)

    # Time loop
     @time for it=1:Nt

        # Compute Ricker function
        t          += Δt
        a           = Ricker(t, t₀, 𝑓₀)
        f_ext[isrc] = ρ[isrc]*a/Δx/Δx

        # Velocity gradient components
        @. ∂Vx∂x[2:end-1] = (V.x[2:end] - V.x[1:end-1])/Δx
        
        # Divergence
        @. ∇V   = ∂Vx∂x

        # Deviatoric strain rate 
        @. ε̇.xx = ∂Vx∂x - 1/3*∇V
      
        # Stress update
        @. τ.xx = f_shear(G,η,Δt)*Δt*(ε̇.xx) + f_relax(G,η,Δt)*τ.xx

        # Pressure update 
        @. P    = P - Δt*f_bulk(K)*∇V

        # Linear momentum balance
        @. V.x[2:end-1] = V.x[2:end-1] + Δt/ρ[2:end-1]*((τ.xx[3:end-1]-τ.xx[2:end-2])/Δx - (P[3:end-1]-P[2:end-2])/Δx - f_ext[2:end-1])

        # Absorbing boundary Cerjean et al. (1985)
        @.  V.x  = V.x  * bc_filtW_v 
        @.  P    = P    * bc_filtW_c 
        @.  τ.xx = τ.xx * bc_filtW_c 
        @.  V.x  = V.x  * bc_filtE_v 
        @.  P    = P    * bc_filtE_c 
        @.  τ.xx = τ.xx * bc_filtE_c 

        # Visualisation
        if mod(it, Nout)==0 && visu==true
            display(plot(xv, V.x, ylim=(-2e-3, 2e-3)))
            #display(plot(xv, V.x))
            sleep(0.1)
        end
    end
    #display(plot(xv, V.x, ylim=(-2e-3, 2e-3)))
    
end



MainSource()
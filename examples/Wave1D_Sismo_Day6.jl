using SeismicQ, Plots, FastBroadcast 

function MainSource()
    
    # Spatial extent
    Lx   = 50.0

    # Discretization
    Ncx = 200
    Δx  = Lx/Ncx
    xv  = LinRange(0,Lx,Ncx+1)
    xc  = LinRange(0-Δx/2,Lx+Δx/2,Ncx+2)

    # Source parameters
    𝑓₀   = 200     # Central frequency of the source [Hz]
    t₀   = 1.2/𝑓₀
    σ₀   = Lx/100
    isrc = Int((Ncx/2)+1)
    x₀   = xv[isrc]

    # Mechanical parameters 
    ρ₀      = 1500.0
    K₀      = 1.e9
    G₀      = 1.e1
    De_s    = 1e-1 # Shear Deborah number
    ηₘ₀     = De_s*G₀ / 𝑓₀
    Fb_b    = 2e-1 # Bulk Fatboy number
    ηₖ₀     = Fb_b*K₀ / 𝑓₀
    DevRheo = :MaxwellVE #:Elastic or :MaxwellVE
    VolRheo = :KelvinVE  #:Elastic or :KelvinVE   
    #DevRheo = :Elastic #or :MaxwellVE
    #VolRheo = :Elastic #or :KelvinVE  


    # Time domain
    c_eff = sqrt((K₀*(1+Fb_b)+4/3*G₀)/ρ₀) 
    Δt    = min(1e10, Δx/c_eff/5) # Courant criteria from wavespeed
    Nt    = 160
    Nout  = 1
    t     = 0.0#-t₀
   
    # Parameters for Sismo.
    Xs   = 25:0.5:50   # x_coordinates [m]
    Ns   = size(Xs,1)
    ds   = zeros(size(Xs))
    @. ds   = abs(Xs-xv[isrc])
    velocity_matrix = zeros(Ns, Nt)
    time = zeros(Nt)
    
    # Storage on centers # +2 for ghost nodes for BCs
    szv   = (Ncx+1,)
    szc   = (Ncx+2,)
    # Storage on centroids 
    K     = ones(szc)*K₀ 
    G     = ones(szc)*G₀
    ηₘ    = ones(szc)*ηₘ₀
    ηₖ    = ones(szc)*ηₖ₀
    ε̇     = ( xx=zeros(szc), yy=zeros(szc), zz=zeros(szc), xy=zeros(szc), yz=zeros(szc), xz=zeros(szc) )  
    ∇V    = zeros(Ncx+2)
    ∇V0   = zeros(Ncx+2)
    P     = zeros(Ncx+2)
    P0    = zeros(Ncx+2)
    τ     = ( xx=zeros(szc), yy=zeros(szc), zz=zeros(szc), xy=zeros(szc), yz=zeros(szc), xz=zeros(szc) )  
    τ0    = ( xx=zeros(szc), yy=zeros(szc), zz=zeros(szc), xy=zeros(szc), yz=zeros(szc), xz=zeros(szc) )  
    ∂Vx∂x = zeros(szc)
    # Storage on vertices
    V     = ( x=zeros(szv), y=zeros(szv), z=zeros(szv))
    ρ     = ones(szv)*ρ₀ 
    f_ext = zeros(szv)

    # BC
    Lbc        = 2.
    bc_filtW_v = 1.0 .- exp.(-(xv.-0Lx).^2/Lbc.^2)
    bc_filtW_c = 1.0 .- exp.(-(xc.-0Lx).^2/Lbc.^2)
    bc_filtE_v = 1.0 .- exp.(-(xv.- Lx).^2/Lbc.^2)
    bc_filtE_c = 1.0 .- exp.(-(xc.- Lx).^2/Lbc.^2)

    # Select deviatoric rheology
    if DevRheo == :Elastic
        dev = (G,Δt)
     elseif DevRheo == :MaxwellVE
        dev = (G,ηₘ,Δt)
    end

    # Select volumetric rheology
    if VolRheo == :Elastic
        vol = (K,Δt)
    elseif VolRheo == :KelvinVE
        vol = (K,ηₖ,Δt)
    end


    # Time loop
    @time for it=1:Nt

        P0    .= P
        @. P0     = P0+χb(dev...)
        τ0.xx .= τ.xx
        @. τ0.xx = τ0.xx-χs(dev...)

        # Compute Ricker function
        t          += Δt
        @. f_ext    = ρ.*Ricker.(xv, x₀, t, t₀, 𝑓₀, σ₀)

        # Velocity gradient components
        @. ∂Vx∂x[2:end-1] = (V.x[2:end] - V.x[1:end-1])/Δx
        
        # Divergence
        @. ∇V   = ∂Vx∂x

        # Deviatoric strain rate 
        @. ε̇.xx = ∂Vx∂x - 1/3*∇V
      
        # Deviatoric stress update
        
        @. τ.xx = ηs(dev...)*ε̇.xx + τ0.xx *θs(dev...)        

        # Pressure update 
        
        @. P    = -ηb(vol...)*∇V + θb(vol...)*P0
        
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
        if mod(it, Nout)==0
            p = plot(xv, V.x,ylim=(-2e-4, 2e-4))
            display(p)
            sleep(0.1)
        end

        # Extract sismo data:
        time[it] = t
        @. velocity_matrix[:,it] = V.x[Int(Xs[:]/Δx)+1]
        
    end
    
    # Visualization Receiver Gather:
    valim = max(abs(maximum(velocity_matrix)),abs(minimum(velocity_matrix)))
    p = heatmap(ds,time,velocity_matrix',color=palette(:RdBu,100,rev=true),
    title="Receiver gather", xlabel="distance to the source [m]",
    ylabel="time [s]",yflip=true,clim=(-valim,+valim))
    display(p)

end

MainSource()
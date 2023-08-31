using SeismicQ, FastBroadcast, GLMakie, Printf, Colors, ColorSchemes, MathTeXEngine
Makie.update_theme!(fonts = (regular = texfont(), bold = texfont(:bold), italic = texfont(:italic)))

function MainSource()
    visu     = true
    printfig = true  # print figures to disk
    path     = "./runs/"
    juliadivcmap    = zeros(RGB{Float64}, 5)
    juliadivcmap[1] = RGBA{Float64}(0/255,150/255,0/255, 1.)  
    juliadivcmap[2] = RGBA{Float64}(0/255,0/255,200/255, 1.)  
    juliadivcmap[3] = RGBA{Float64}(255/255,255/255,255/255, 1.) 
    juliadivcmap[4] = RGBA{Float64}(150/255,0/255,150/255, 1.) 
    juliadivcmap[5] = RGBA{Float64}(200/255,0/255,0/255, 1.)
    wave_colors     = cgrad(juliadivcmap, length(juliadivcmap), categorical=false, rev=false)
    
    # Spatial extent
    l  = (x = 100, y = 25)

    # Discretization
    Nc  = (x = 200, y = 200) 
    Δ   = (x = l.x/Nc.x, y = l.y/Nc.y, z=1.0)
    X   = (v = (x= LinRange(0,l.x,Nc.x+1)             , y= LinRange(0,l.y,Nc.y+1)),
           c = (x= LinRange(0-Δ.x/2,l.x+Δ.x/2,Nc.x+2) , y= LinRange(0-Δ.y/2,l.y+Δ.y/2,Nc.y+2)),
           i = (x= LinRange(0,l.x,Nc.x+1)             , y= LinRange(0-Δ.y/2,l.y+Δ.y/2,Nc.y+2)),
           j = (x= LinRange(0-Δ.x/2,l.x+Δ.x/2,Nc.x+2) , y= LinRange(0,l.y,Nc.y+1))) 
        
    # Source parameters
    𝑓₀   = 100   # Central frequency of the source [Hz]
    t₀   = 1.2/𝑓₀
    σ₀   = l.x/100
    x₀   = l.x/2
    y₀   = l.y/2
    src  = (i=Int((Nc.x/2)+1),j=Int((Nc.y/2)+1))
    facS = (v=(x=0.0,y=1.0,z=1.0),c=(x=0.0,y=1.0,z=1.0))
    
  


    # Mechanical parameters 
    ρ₀      = 1500.0
    K₀      = 1.e9
    G₀      = 1.e8
    De_s    = 1e-2 # Shear Deborah number
    ηₘ₀     = De_s*G₀ / 𝑓₀
    Fb_b    = 1e-2 # Bulk Fatboy number
    ηₖ₀     = Fb_b*K₀ / 𝑓₀
   # DevRheo = :MaxwellVE #:Elastic or :MaxwellVE
   # VolRheo = :KelvinVE  #:Elastic or :KelvinVE 
    
    DevRheo = :Elastic #or :MaxwellVE
    VolRheo = :Elastic #or :KelvinVE 

    # Time domain
    c_eff = sqrt((K₀*(1+Fb_b)+4/3*G₀)/ρ₀) 
    Δt    = min(1e10, 0.1*Δ.x/c_eff, 0.1*Δ.y/c_eff ) # Courant criteria from wavespeed
    Nt    = 10001
    Nout  = 250
    t     = -t₀


    # Parameters for Sismo.
    Noutsismo       = 10
    Ntsismo         = Int((Nt-1)/Noutsismo)
    Xs              = LinRange(0,l.x,Nc.x+1)*ones(1,Ntsismo)   # x_coordinates [m]
    Ns              = size(Xs,1)
    velocity_matrix = (x = zeros(Ns, Ntsismo), y=zeros(Ns, Ntsismo), z=zeros(Ns, Ntsismo)) 
    arrival_time    = zeros(Ns,Ntsismo)
    ksismo          = 0
   
    # Storage on centers # +2 for ghost nodes for BCs
    szv   = (Nc.x+1, Nc.y+1)
    szc   = (Nc.x+2, Nc.y+2)
    szi   = (Nc.x+1, Nc.y+2)
    szj   = (Nc.x+2, Nc.y+1)
    # Storage on i and j meshes
    K     = (i= ones(szi)*K₀,  j= ones(szj)*K₀ ) 
    G     = (i= ones(szi)*G₀,  j= ones(szj)*G₀ ) 
    ηₘ    = (i= ones(szi)*ηₘ₀ , j= ones(szj)*ηₘ₀)
    ηₖ    = (i= ones(szi)*ηₖ₀ , j= ones(szj)*ηₖ₀ )
    ∇V    = (i=zeros(szi),     j=zeros(szj))
    P     = (i=zeros(szi),     j=zeros(szj))
    P0   = (i=zeros(szi),     j=zeros(szj))
    L     = (i=(xx=zeros(szi), xy=zeros(szi), yx=zeros(szi), yy=zeros(szi),zx=zeros(szi),zy=zeros(szi)),
             j=(xx=zeros(szj), xy=zeros(szj), yx=zeros(szj), yy=zeros(szj),zx=zeros(szj),zy=zeros(szj)))
    ε̇     = (i=(xx=zeros(szi), yy=zeros(szi), zz=zeros(szi), xy=zeros(szi),xz=zeros(szi),yz=zeros(szi)),
             j=(xx=zeros(szj), yy=zeros(szj), zz=zeros(szj), xy=zeros(szj),xz=zeros(szj),yz=zeros(szj))) 
    τ     = (i=(xx=zeros(szi), yy=zeros(szi), zz=zeros(szi), xy=zeros(szi),xz=zeros(szi),yz=zeros(szi)),
             j=(xx=zeros(szj), yy=zeros(szj), zz=zeros(szj), xy=zeros(szj),xz=zeros(szj),yz=zeros(szj))) 
    τ0    = (i=(xx=zeros(szi), yy=zeros(szi), zz=zeros(szi), xy=zeros(szi),xz=zeros(szi),yz=zeros(szi)),
             j=(xx=zeros(szj), yy=zeros(szj), zz=zeros(szj), xy=zeros(szj),xz=zeros(szj),yz=zeros(szj)))                   

    # Storage on v and c meshes
    V     = ( v=(x=zeros(szv), y=zeros(szv), z=zeros(szv)),
              c=(x=zeros(szc), y=zeros(szc), z=zeros(szc)))

    ρ     = (v=ones(szv)*ρ₀, c=ones(szc)*ρ₀)
    f_ext = (v=zeros(szv)  , c=zeros(szc))
    # BC
    Lbc        = 1.
    # BC on v and c mesh
    bc_filt_V   = (v=Cerjean2D(X.v,Lbc,l,Δ),c=Cerjean2D(X.c,Lbc,l,Δ))
    bc_filt_tau = (i=Cerjean2D(X.i,Lbc,l,Δ),j=Cerjean2D(X.j,Lbc,l,Δ))

    # Compute Ricker function with 2D spatial support
    f_ext = (v=zeros(szv)  , c=zeros(szc))
    xc2d   = X.c.x * ones(size( X.c.y))'
    yc2d   = ones(size( X.c.x)) * X.c.y'
    xv2d   = X.v.x * ones(size( X.v.y))'
    yv2d   = ones(size( X.v.x)) * X.v.y'

    # Select deviatoric rheology
    if DevRheo == :Elastic
        devi = (G.i,Δt)
        devj = (G.j,Δt)
    elseif DevRheo == :MaxwellVE
        devi = (G.i,ηₘ.i,Δt)
        devj = (G.j,ηₘ.j,Δt)
    end

    # Select volumetric rheology
    if VolRheo == :Elastic
        voli = (K.i,Δt)
        volj = (K.j,Δt)
    elseif VolRheo == :KelvinVE
        voli = (K.i,ηₖ.i,Δt)
        volj = (K.j,ηₖ.j,Δt)
    end

    # Time loop
    @views @time for it=1:Nt

        # Update Time
        t += Δt
        P0.i .= P.i
        P0.j .= P.j

        τ0.i.xx .= τ.i.xx
        τ0.i.xy .= τ.i.xy
        τ0.i.xz .= τ.i.xz
        τ0.i.yy .= τ.i.yy
        τ0.i.zz .= τ.i.zz
        τ0.i.yz .= τ.i.yz

        τ0.j.xx .= τ.j.xx
        τ0.j.xy .= τ.j.xy
        τ0.j.xz .= τ.j.xz
        τ0.j.yy .= τ.j.yy
        τ0.j.zz .= τ.j.zz
        τ0.j.yz .= τ.j.yz 
        
        # 2D Ricker with spatial support
        @.. f_ext.c = ρ.c*Ricker.( xc2d, x₀, yc2d, y₀, t, t₀, 𝑓₀, σ₀)
        @.. f_ext.v = ρ.v*Ricker.( xv2d, x₀, yv2d, y₀, t, t₀, 𝑓₀, σ₀)
        
        # Inherited pressure (remove the instantaneous viscous contribution )
       
        @.. P0.i = P0.i + χb(voli...)*∇V.i 
        @.. P0.j = P0.j + χb(volj...)*∇V.j

        # Inherited deviatoric stress (remove the instantaneous viscous contribution)
        @.. τ0.i.xx= τ0.i.xx - χs(devi...)*ε̇.i.xx
        @.. τ0.i.xy= τ0.i.xy - χs(devi...)*ε̇.i.xy
        @.. τ0.i.xz= τ0.i.xz - χs(devi...)*ε̇.i.xz
        @.. τ0.i.yy= τ0.i.yy - χs(devi...)*ε̇.i.yy
        @.. τ0.i.zz= τ0.i.zz - χs(devi...)*ε̇.i.zz
        @.. τ0.i.yz= τ0.i.yz - χs(devi...)*ε̇.i.yz
       
        @.. τ0.j.xx= τ0.j.xx - χs(devj...)*ε̇.j.xx
        @.. τ0.j.xy= τ0.j.xy - χs(devj...)*ε̇.j.xy
        @.. τ0.j.xz= τ0.j.xz - χs(devj...)*ε̇.j.xz
        @.. τ0.j.yy= τ0.j.yy - χs(devj...)*ε̇.j.yy
        @.. τ0.j.zz= τ0.j.zz - χs(devj...)*ε̇.j.zz
        @.. τ0.j.yz= τ0.j.yz - χs(devj...)*ε̇.j.yz

        # Velocity gradient components
        @.. L.i.xx[:,2:end-1] = (V.c.x[2:end,2:end-1] - V.c.x[1:end-1,2:end-1])/Δ.x
        @.. L.j.xx[2:end-1,:] = (V.v.x[2:end,:] - V.v.x[1:end-1,:])/Δ.x

        @.. L.i.yx[:,2:end-1] = (V.c.y[2:end,2:end-1] - V.c.y[1:end-1,2:end-1])/Δ.x
        @.. L.j.yx[2:end-1,:] = (V.v.y[2:end,:] - V.v.y[1:end-1,:])/Δ.x

        @.. L.i.yy[:,2:end-1] = (V.v.y[:,2:end] - V.v.y[:,1:end-1])/Δ.y
        @.. L.j.yy[2:end-1,:] = (V.c.y[2:end-1,2:end] - V.c.y[2:end-1,1:end-1])/Δ.y

        @.. L.i.xy[:,2:end-1] = (V.v.x[:,2:end] - V.v.x[:,1:end-1])/Δ.y
        @.. L.j.xy[2:end-1,:] = (V.c.x[2:end-1,2:end] - V.c.x[2:end-1,1:end-1])/Δ.y

        @.. L.i.zy[:,2:end-1] = (V.v.z[:,2:end] - V.v.z[:,1:end-1])/Δ.y
        @.. L.j.zy[2:end-1,:] = (V.c.z[2:end-1,2:end] - V.c.z[2:end-1,1:end-1])/Δ.y

        @.. L.i.zx[:,2:end-1] = (V.c.z[2:end,2:end-1] - V.c.z[1:end-1,2:end-1])/Δ.x
        @.. L.j.zx[2:end-1,:] = (V.v.z[2:end,:] - V.v.z[1:end-1,:])/Δ.x
        
        # Divergence
        @.. ∇V.i   = L.i.xx + L.i.yy
        @.. ∇V.j   = L.j.xx + L.j.yy

        # Deviatoric strain rate 
        @.. ε̇.i.xx = L.i.xx - 1//3*∇V.i
        @.. ε̇.j.xx = L.j.xx - 1//3*∇V.j

        @.. ε̇.i.yy = L.i.yy - 1//3*∇V.i
        @.. ε̇.j.yy = L.j.yy - 1//3*∇V.j

        @.. ε̇.i.zz = - 1//3*∇V.i
        @.. ε̇.j.zz = - 1//3*∇V.j

        @.. ε̇.i.xy = 1//2*(L.i.xy + L.i.yx)
        @.. ε̇.j.xy = 1//2*(L.j.xy + L.j.yx)
        
        # in 2D Lxz and Lyz are zero 
        @.. ε̇.i.xz = 1//2*(L.i.zx)
        @.. ε̇.j.xz = 1//2*(L.j.zx)

        @.. ε̇.i.yz = 1//2*(L.i.zy)
        @.. ε̇.j.yz = 1//2*(L.j.zy)
      
        # Stress update
        @.. τ.i.xx = ηs(devi...)*(ε̇.i.xx) + θs(devi...)*τ0.i.xx
        @.. τ.j.xx = ηs(devj...)*(ε̇.j.xx) + θs(devj...)*τ0.j.xx

        @.. τ.i.yy = ηs(devi...)*(ε̇.i.yy) + θs(devi...)*τ0.i.yy
        @.. τ.j.yy = ηs(devj...)*(ε̇.j.yy) + θs(devj...)*τ0.j.yy
        
        @.. τ.i.zz = ηs(devi...)*(ε̇.i.zz) + θs(devi...)*τ0.i.zz
        @.. τ.j.zz = ηs(devj...)*(ε̇.j.zz) + θs(devj...)*τ0.j.zz
        
        @.. τ.i.xy = ηs(devi...)*(ε̇.i.xy) + θs(devi...)*τ0.i.xy
        @.. τ.j.xy = ηs(devj...)*(ε̇.j.xy) + θs(devj...)*τ0.j.xy
        
        @.. τ.i.xz = ηs(devi...)*(ε̇.i.xz) + θs(devi...)*τ0.i.xz
        @.. τ.j.xz = ηs(devj...)*(ε̇.j.xz) + θs(devj...)*τ0.j.xz
        
        @.. τ.i.yz = ηs(devi...)*(ε̇.i.yz) + θs(devi...)*τ0.i.yz
        @.. τ.j.yz = ηs(devj...)*(ε̇.j.yz) + θs(devj...)*τ0.j.yz

       

        # Pressure update 

        @.. P.i    = θb(voli...)*P0.i - ηb(voli...)*∇V.i 
        @.. P.j    = θb(volj...)*P0.j - ηb(volj...)*∇V.j 


        τ.j.xy[:,end-1:end] .= 0.
        τ.j.yz[:,end-1:end] .= 0.
        τ.j.yy[:,end-1:end] .= P.j[:,end-1:end]

        # Linear momentum balance
        @.. V.v.x[2:end-1,2:end-1] = (V.v.x[2:end-1,2:end-1] 
                                    + Δt/ρ.v[2:end-1,2:end-1]
                                    *((τ.j.xx[3:end-1,2:end-1]-τ.j.xx[2:end-2,2:end-1])/Δ.x
                                    + (τ.i.xy[2:end-1,3:end-1]-τ.i.xy[2:end-1,2:end-2])/Δ.y 
                                    - (P.j[3:end-1,2:end-1]-P.j[2:end-2,2:end-1])/Δ.x 
                                    - facS.v.x*f_ext.v[2:end-1,2:end-1]))
        @.. V.c.x[2:end-1,2:end-1] = (V.c.x[2:end-1,2:end-1] 
                                    + Δt/ρ.c[2:end-1,2:end-1]
                                    *((τ.i.xx[2:end,2:end-1]-τ.i.xx[1:end-1,2:end-1])/Δ.x
                                    + (τ.j.xy[2:end-1,2:end]-τ.j.xy[2:end-1,1:end-1])/Δ.y
                                    - (P.i[2:end,2:end-1]-P.i[1:end-1,2:end-1])/Δ.x 
                                    - facS.c.x*f_ext.c[2:end-1,2:end-1]))                            

        @.. V.v.y[2:end-1,2:end-1] = (V.v.y[2:end-1,2:end-1] 
                                    + Δt/ρ.v[2:end-1,2:end-1]
                                    *((τ.j.xy[3:end-1,2:end-1]-τ.j.xy[2:end-2,2:end-1])/Δ.x
                                    + (τ.i.yy[2:end-1,3:end-1]-τ.i.yy[2:end-1,2:end-2])/Δ.y 
                                    - (P.i[2:end-1,3:end-1]-P.i[2:end-1,2:end-2])/Δ.y 
                                    - facS.v.y*f_ext.v[2:end-1,2:end-1]))
        
        @.. V.c.y[2:end-1,2:end-1] = (V.c.y[2:end-1,2:end-1] 
                                    + Δt/ρ.c[2:end-1,2:end-1]
                                    *((τ.i.xy[2:end,2:end-1]-τ.i.xy[1:end-1,2:end-1])/Δ.x
                                    + (τ.j.yy[2:end-1,2:end]-τ.j.yy[2:end-1,1:end-1])/Δ.y 
                                    - (P.j[2:end-1,2:end]-P.j[2:end-1,1:end-1])/Δ.y 
                                    - facS.c.y*f_ext.c[2:end-1,2:end-1]))   

        # the two terms in dPdz and dtauzzdz  cancel in linear elastic case ... but i am not sure with other rheologies so I have left them 
        @.. V.v.z[2:end-1,2:end-1] = (V.v.z[2:end-1,2:end-1] 
                                    + Δt/ρ.v[2:end-1,2:end-1]
                                    *((τ.j.xz[3:end-1,2:end-1]-τ.j.xz[2:end-2,2:end-1])/Δ.x
                                    + (τ.i.yz[2:end-1,3:end-1]-τ.i.yz[2:end-1,2:end-2])/Δ.y 
                                    - facS.v.z* f_ext.v[2:end-1,2:end-1]))
        
        @.. V.c.z[2:end-1,2:end-1] = (V.c.z[2:end-1,2:end-1] 
                                    + Δt/ρ.c[2:end-1,2:end-1]
                                    *((τ.i.xz[2:end,2:end-1]-τ.i.xz[1:end-1,2:end-1])/Δ.x
                                    + (τ.j.yz[2:end-1,2:end]-τ.j.yz[2:end-1,1:end-1])/Δ.y 
                                    - facS.c.z*f_ext.c[2:end-1,2:end-1]))  
                                    
        # free surface BC at top 
        @.. V.v.y[2:end-1,end] = V.v.y[2:end-1,end-1]
        @.. V.v.x[2:end-1,end] = V.v.x[2:end-1,end-1]
        @.. V.v.z[2:end-1,end] = V.v.z[2:end-1,end-1]

    
        # Absorbing boundary Cerjean et al. (1985)
        @..  V.v.x  = V.v.x  * bc_filt_V.v 
        @..  V.v.y  = V.v.y  * bc_filt_V.v 
        @..  V.v.z  = V.v.z  * bc_filt_V.v
        @..  V.c.x  = V.c.x  * bc_filt_V.c 
        @..  V.c.y  = V.c.y  * bc_filt_V.c 
        @..  V.c.z  = V.c.z  * bc_filt_V.c  

        @..  P.i    = P.i    *  bc_filt_tau.i 
        @..  τ.i.xx = τ.i.xx *  bc_filt_tau.i
        @..  τ.i.yy = τ.i.yy *  bc_filt_tau.i
        @..  τ.i.zz = τ.i.zz *  bc_filt_tau.i
        @..  τ.i.xy = τ.i.xy *  bc_filt_tau.i
        @..  τ.i.xz = τ.i.xz *  bc_filt_tau.i
        @..  τ.i.yz = τ.i.yz *  bc_filt_tau.i

        @..  P.j    = P.j    *  bc_filt_tau.j 
        @..  τ.j.xx = τ.j.xx *  bc_filt_tau.j
        @..  τ.j.yy = τ.j.yy *  bc_filt_tau.j
        @..  τ.j.zz = τ.j.zz *  bc_filt_tau.j
        @..  τ.j.xy = τ.j.xy *  bc_filt_tau.j
        @..  τ.j.xz = τ.j.xz *  bc_filt_tau.j
        @..  τ.j.yz = τ.j.yz *  bc_filt_tau.j

        # Visualisation
        if mod(it, Nout)==0 && visu==true

            resol=500 
            f = Figure(resolution = (l.x/l.y*resol*2, resol*2), fontsize=15)

            ax1 = Axis(f[1, 1], aspect=l.x/l.y, title = L" vx on v grid at $t$ = %$(t) [s]", xlabel = L"$x$ [m]", ylabel = L"$y$ [m]")
            hm = GLMakie.heatmap!(ax1, X.v.x, X.v.y, V.v.x, colormap = wave_colors,colorrange=(-3.e-5,3.e-5))
        
            ax2 = Axis(f[2, 1], aspect=l.x/l.y, title = L" vx on c grid at $t$ = %$(t) [s]", xlabel = L"$x$ [m]", ylabel = L"$y$ [m]")
            hm = GLMakie.heatmap!(ax2, X.c.x, X.c.y, V.c.x, colormap = wave_colors,colorrange=(-3.e-5,3.e-5))

            ax3 = Axis(f[1, 2], aspect=l.x/l.y, title = L" vy on v grid at $t$ = %$(t) [s]", xlabel = L"$x$ [m]", ylabel = L"$y$ [m]")
            hm = GLMakie.heatmap!(ax3, X.v.x, X.v.y, V.v.y, colormap = wave_colors,colorrange=(-3.e-5,3.e-5))
        
            ax4 = Axis(f[2, 2], aspect=l.x/l.y, title = L" vy on c grid at $t$ = %$(t) [s]", xlabel = L"$x$ [m]", ylabel = L"$y$ [m]")
            hm = GLMakie.heatmap!(ax4, X.c.x, X.c.y, V.c.y, colormap = wave_colors,colorrange=(-3.e-5,3.e-5))

            ax5 = Axis(f[1, 3], aspect=l.x/l.y, title = L" vz on v grid at $t$ = %$(t) [s]", xlabel = L"$x$ [m]", ylabel = L"$y$ [m]")
            hm = GLMakie.heatmap!(ax5, X.v.x, X.v.y, V.v.z, colormap = wave_colors,colorrange=(-3.e-5,3.e-5))
        
            ax6 = Axis(f[2, 3], aspect=l.x/l.y, title = L" vz on c grid at $t$ = %$(t) [s]", xlabel = L"$x$ [m]", ylabel = L"$y$ [m]")
            hm = GLMakie.heatmap!(ax6, X.c.x, X.c.y, V.c.z, colormap = wave_colors,colorrange=(-3.e-5,3.e-5))

            # colsize!(f.layout, 1, Aspect(1, l.x/l.y))
            GLMakie.Colorbar(f[1, 4], hm, label = "V [m/s]", width = 20, labelsize = 25, ticklabelsize = 14 )
            # GLMakie.colgap!(f.layout, 20)
            display(f)
            sleep(0.1)
            if printfig Print2Disk( f, path, "Vall", it) end
        end
        if mod(it, Noutsismo)==0 && visu==true
         # Extract sismo data:
         ksismo +=1 
         arrival_time[:,ksismo] .= t
         velocity_matrix.x[:,ksismo] .= V.v.x[:,end-1]
         velocity_matrix.y[:,ksismo] .= V.v.y[:,end-1]
         velocity_matrix.z[:,ksismo] .= V.v.z[:,end-1]
        end
    end

    valimx = max(abs(maximum(velocity_matrix.x)),abs(minimum(velocity_matrix.x)))
    valimy = max(abs(maximum(velocity_matrix.y)),abs(minimum(velocity_matrix.y)))
    valimz = max(abs(maximum(velocity_matrix.z)),abs(minimum(velocity_matrix.z)))

     resol=500 
             f = Figure(resolution = (l.x/l.y*resol*3, resol), fontsize=15)

             ax1 = Axis(f[1, 1],  title = L" vx [m/s]", xlabel = L"$x$ [m]", ylabel = L"$t$ [s]")
             hm = GLMakie.heatmap!(ax1,  velocity_matrix.x, colormap = wave_colors,colorrange=(-3.e-5,3.e-5))
        
            ax2 = Axis(f[1, 2],  title = L" vy [m/s]", xlabel = L"$x$ [m]", ylabel = L"$t$ [s]")
             hm = GLMakie.heatmap!(ax2,velocity_matrix.y, colormap = wave_colors,colorrange=(-3.e-5,3.e-5))

            ax3 = Axis(f[1, 3], title = L" vz [m/s]", xlabel = L"$x$ [m]", ylabel = L"$t$ [m]")
            hm = GLMakie.heatmap!(ax3, velocity_matrix.z, colormap = wave_colors,colorrange=(-3.e-5,3.e-5))
        
             GLMakie.Colorbar(f[1, 4], hm, label = "V [m/s]", width = 20, labelsize = 25, ticklabelsize = 14 )
        #  # GLMakie.colgap!(f.layout, 20)
         display(f)

end

function Cerjean2D(X,Lbc,l,Δ)
    return ((1.0 .- exp.(-(X.x*ones(size(X.y))'.-0*l.x).^2/Lbc.^2))
         .*(1.0 .- exp.(-(X.x*ones(size(X.y))' .-  l.x).^2/Lbc.^2))
         .*(1.0 .- exp.(-(ones(size(X.x))*X.y' .-0*l.y).^2/Lbc.^2))
         )
         #.*(1.0 .- exp.(-(ones(size(X.x))*X.y' .-  l.y).^2/Lbc.^2))
end

function Print2Disk( f, path, field, istep; res=4)
    path1 = path*"/_$field/"
    mkpath(path1)
    save(path1*"$field"*@sprintf("%05d", istep)*".png", f, px_per_unit = res) 
end

MainSource()
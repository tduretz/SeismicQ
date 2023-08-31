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
    l  = (x = 100, y = 50)

    # Discretization
    Nc  = (x = 400, y = 100) 
    Δ   = (x = l.x/Nc.x, y = l.y/Nc.y, z=1.0)
    X   = (v = (x= LinRange(0,l.x,Nc.x+1)             , y= LinRange(0,l.y,Nc.y+1)),
           c = (x= LinRange(0-Δ.x/2,l.x+Δ.x/2,Nc.x+2) , y= LinRange(0-Δ.y/2,l.y+Δ.y/2,Nc.y+2)),
           i = (x= LinRange(0,l.x,Nc.x+1)             , y= LinRange(0-Δ.y/2,l.y+Δ.y/2,Nc.y+2)),
           j = (x= LinRange(0-Δ.x/2,l.x+Δ.x/2,Nc.x+2) , y= LinRange(0,l.y,Nc.y+1))) 
        
    # Source parameters
    𝑓₀   = 50  # Central frequency of the source [Hz]
    t₀   = 1.2/𝑓₀
    σ₀   = l.x/100
    x₀   = l.x/2
    y₀   = l.y/2
    src  = (i=Int((Nc.x/2)+1),j=Int((Nc.y/4)+1))
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
    Nout  = 500
    t     = -t₀


    # Parameters for Sismo.
    Noutsismo       = 10
    Ntsismo         = Int((Nt-1)/Noutsismo)
    Xs              = LinRange(0,l.x,Nc.x+1)*ones(1,Ntsismo).*1.0   # x_coordinates [m]
    Ns              = size(Xs,1)
    velocity_matrix = (x = zeros(Ns, Ntsismo).*1.0, y=zeros(Ns, Ntsismo).*1.0, z=zeros(Ns, Ntsismo).*1.0) 
    arrival_time    = zeros(Ns,Ntsismo).*1.0
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
    Lbc        = 3.
    # BC on v and c mesh
    Δ0 = (x=0,y=0)
    Δi = (x=0,y=Δ.y/2)
    Δj = (x=Δ.x/2,y=0)
    bc_filt_V   = (v=Cerjean2D(X.v,Lbc,l,Δ0),c=Cerjean2D(X.c,Lbc,l,Δ))
    bc_filt_tau = (i=Cerjean2D(X.i,Lbc,l,Δi),j=Cerjean2D(X.j,Lbc,l,Δj))

    # Compute Ricker function with 2D spatial support
    f_ext = (v=zeros(szv)  , c=zeros(szc))
    xc2d   = X.c.x * ones(size( X.c.y))'
    yc2d   = ones(size( X.c.x)) * X.c.y'
    xv2d   = X.v.x * ones(size( X.v.y))'
    yv2d   = ones(size( X.v.x)) * X.v.y'

    # Select deviatoric rheology
    if DevRheo == :Elastic
        dev = (i=(G.i,Δt), j=(G.j,Δt)) 
    elseif DevRheo == :MaxwellVE
        dev = (i=(G.i,ηₘ.i,Δt), j= (G.j,ηₘ.j,Δt))
    end

    # Select volumetric rheology
    if VolRheo == :Elastic
        vol = (i=(K.i,Δt),j=(K.j,Δt))
    elseif VolRheo == :KelvinVE
        vol = (i=(K.i,ηₖ.i,Δt),j=(K.j,ηₖ.j,Δt))
    end


    # Time loop
    @views @time for it=1:Nt

        # Update Time
        t += Δt
        P0.i .= P.i
        P0.j .= P.j
        for grid=1:2, comp=1:length(τ0.i) 
            τ0[grid][comp] .= τ[grid][comp]
        end
    
        
        # 2D Ricker with spatial support
        @.. f_ext.c = ρ.c*Ricker.( xc2d, x₀, yc2d, y₀, t, t₀, 𝑓₀, σ₀)
        @.. f_ext.v = ρ.v*Ricker.( xv2d, x₀, yv2d, y₀, t, t₀, 𝑓₀, σ₀)
        
        # Inherited pressure (remove the instantaneous viscous contribution )
       
        @.. P0.i = P0.i + χb(vol.i...)*∇V.i 
        @.. P0.j = P0.j + χb(vol.j...)*∇V.j

        # Inherited deviatoric stress (remove the instantaneous viscous contribution)
        for grid=1:2, comp=1:length(τ0.i) 
            @.. τ0[grid][comp]= τ0[grid][comp] - χs(dev[grid]...)*ε̇[grid][comp]
        end
        

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
  
        
        # L.j.xy[:,end] .= (-L.j.yx[:,end] .* ηs(dev.j...)[:,end] .- theta_s[:,end]  .* τ0.j.xy[:,end]) ./ ηs(dev.j...)[:,end]
        # L.j.zy[:,end] .= -θs(dev.j...)[:,end].*  τ0.j.yz[:,end] ./ ηs(dev.j...)[:,end]
        L.j.yy[:,end] .= (- 3.0 * L.j.xx[:,end] .* ηb(vol.j...)[:,end] + 2.0 * L.j.xx[:,end] .* ηs(dev.j...)[:,end] 
                         + 3.0 * P0.j[:,end]   .* θb(vol.j...)[:,end] - 3.0 * θs(dev.j...)[:,end]   .* τ0.j.yy[:,end]) ./
                         (3.0 * ηb(vol.j...)[:,end] + 4.0 * ηs(dev.j...)[:,end])

        L.i.xy[:,end-1] .= (-L.i.yx[:,end-1] .* ηs(dev.i...)[:,end-1] .- θs(dev.i...)[:,end-1].*τ0.i.xy[:,end-1]) ./ ηs(dev.i...)[:,end-1]
        L.i.zy[:,end-1] .= -θs(dev.i...)[:,end-1].*  τ0.i.yz[:,end-1] ./ ηs(dev.i...)[:,end-1]
        L.i.yy[:,end-1] .= (- 3.0 * L.i.xx[:,end-1] .* ηb(vol.i...)[:,end-1] + 2.0 * L.i.xx[:,end-1] .* ηs(dev.i...)[:,end-1]
                         + 3.0 * P0.i[:,end-1]   .* θb(vol.i...)[:,end-1]  - 3.0 * θs(dev.i...)[:,end-1]  .* τ0.i.yy[:,end-1]) ./
                         (3.0 * ηb(vol.i...)[:,end-1]  + 4.0 * ηs(dev.i...)[:,end-1])
        
        
        
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

        for grid=1:2, comp=1:length(τ0.i) 
             @.. τ[grid][comp] = ηs(dev[grid]...)*(ε̇[grid][comp]) + θs(dev[grid]...)*τ0[grid][comp]
        end
       

        
        # Pressure update 
        for grid=1:2
            @.. P[grid]  = θb(vol[grid]...)*P0[grid] - ηb(vol[grid]...)*∇V[grid]
        end 
        #@.. P.j    = θb(vol.j...)*P0.j - ηb(vol.j...)*∇V.j 


        τ.j.xy[:,end] .= 0.
        τ.j.yz[:,end] .= 0.
        τ.j.yy[:,end] .= P.j[:,end]

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
        # @.. V.v.y[2:end-1,end] = V.v.y[2:end-1,end-1]
        # @.. V.v.x[2:end-1,end] = V.v.x[2:end-1,end-1]
        # @.. V.v.z[2:end-1,end] = V.v.z[2:end-1,end-1]

       
        @.. V.v.x[2:end-1,end] = (V.v.x[2:end-1,end] 
                                 + Δt/ρ.v[2:end-1,end]
                                 *((τ.j.xx[3:end-1,end]-τ.j.xx[2:end-2,end])/Δ.x
                                 + (τ.i.xy[2:end-1,end]-τ.i.xy[2:end-1,end-1])/Δ.y 
                                 - (P.j[3:end-1,end]-P.j[2:end-2,end])/Δ.x 
                                 - facS.v.x*f_ext.v[2:end-1,end]))

        @.. V.v.y[2:end-1,end] = (V.v.y[2:end-1,end] 
                                 + Δt/ρ.v[2:end-1,end]
                                 *((τ.j.xy[3:end-1,end]-τ.j.xy[2:end-2,end])/Δ.x
                                 + (τ.i.yy[2:end-1,end]-τ.i.yy[2:end-2,end-1])/Δ.y 
                                 - (P.i[2:end-1,end]   -P.i[2:end-1,end-1])/Δ.y 
                                 - facS.v.y*f_ext.v[2:end-1,end]))

        @.. V.v.z[2:end-1,end] = (V.v.z[2:end-1,end] 
                                    + Δt/ρ.v[2:end-1,end]
                                    *((τ.j.xz[3:end-1,end]-τ.j.xz[2:end-2,end])/Δ.x
                                    + (τ.i.yz[2:end-1,end]-τ.i.yz[2:end-1,end-1])/Δ.y 
                                    - facS.v.z* f_ext.v[2:end-1,end]))


    
        # Absorbing boundary Cerjean et al. (1985)
        for grid=1:2, comp=1:length(V.v) 
            @..  V[grid][comp] = V[grid][comp]  * bc_filt_V[grid] 
        end 
 

        @..  P.i    = P.i    *  bc_filt_tau.i 
        @..  P.j    = P.j    *  bc_filt_tau.j 
        for grid=1:2, comp=1:length(τ.i) 
            @..  τ[grid][comp] *= bc_filt_V[grid] 
        end 


        # Visualisation
        if mod(it, Nout)==0 && visu==true

            resol=200 
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
             f = Figure(resolution = (resol,resol), fontsize=15)

             ax1 = Axis(f[1, 1],  title = L" vx [m/s]", xlabel = L"$x$ [m]", ylabel = L"$t$ [s]")
             hm = GLMakie.heatmap!(ax1, Xs[:,1], -arrival_time[1,:], velocity_matrix.x, colormap = wave_colors,colorrange=(-valimx,valimx))
             GLMakie.Colorbar(f[1, 2], hm, label = "V [m/s]", width = 20, labelsize = 25, ticklabelsize = 14 )
           
            ax2 = Axis(f[2, 1],  title = L" vy [m/s]", xlabel = L"$x$ [m]", ylabel = L"$t$ [s]")
             hm = GLMakie.heatmap!(ax2,Xs[:,1], -arrival_time[1,:],velocity_matrix.y, colormap = wave_colors,colorrange=(-valimy,valimy))
             GLMakie.Colorbar(f[2, 2], hm, label = "V [m/s]", width = 20, labelsize = 25, ticklabelsize = 14 )
           
             ax3 = Axis(f[3, 1], title = L" vz [m/s]", xlabel = L"$x$ [m]", ylabel = L"$t$ [m]")
            hm = GLMakie.heatmap!(ax3, Xs[:,1], -arrival_time[1,:],  velocity_matrix.z, colormap = wave_colors,colorrange=(-valimz,valimz))
             GLMakie.Colorbar(f[3, 2], hm, label = "V [m/s]", width = 20, labelsize = 25, ticklabelsize = 14 )
        #  # GLMakie.colgap!(f.layout, 20)
         display(f)

end

function Cerjean2D(X,Lbc,l,Δ)
    return ((1.0 .- exp.(-(X.x*ones(size(X.y))'.- (0*l.x+Δ.x/2)).^2/Lbc.^2))
         .*(1.0 .- exp.(-(X.x*ones(size(X.y))' .- (l.x-Δ.x/2)).^2/Lbc.^2))
         .*(1.0 .- exp.(-(ones(size(X.x))*X.y' .- (0*l.y+Δ.y/2)).^2/Lbc.^2))
         )
         #.*(1.0 .- exp.(-(ones(size(X.x))*X.y' .-  l.y).^2/Lbc.^2))
end

function Print2Disk( f, path, field, istep; res=4)
    path1 = path*"/_$field/"
    mkpath(path1)
    save(path1*"$field"*@sprintf("%05d", istep)*".png", f, px_per_unit = res) 
end

MainSource()
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
    X   = (x= ( v = LinRange(0,l.x,Nc.x+1)            ,c= LinRange(0-Δ.x/2,l.x+Δ.x/2,Nc.x+2),
                i = LinRange(0,l.x,Nc.x+1)            ,j= LinRange(0-Δ.x/2,l.x+Δ.x/2,Nc.x+2)),
           y= ( v = LinRange(0,l.y,Nc.y+1)            ,c= LinRange(0-Δ.y/2,l.y+Δ.y/2,Nc.y+2),
                i = LinRange(0-Δ.y/2,l.y+Δ.y/2,Nc.y+2),j= LinRange(0,l.y,Nc.y+1)))         
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
    topbc=(i=szi[2]-1,j=szj[2])
    # Storage on i and j meshes
    K     = (i= ones(szi)*K₀,  j= ones(szj)*K₀ ) 
    G     = (i= ones(szi)*G₀,  j= ones(szj)*G₀ ) 
    ηₘ    = (i= ones(szi)*ηₘ₀ , j= ones(szj)*ηₘ₀)
    ηₖ    = (i= ones(szi)*ηₖ₀ , j= ones(szj)*ηₖ₀ )
    ∇V    = (i=zeros(szi),     j=zeros(szj))
    P     = (i=zeros(szi),     j=zeros(szj))
    P0   = (i=zeros(szi),     j=zeros(szj))
    # L     = (i=(xx=zeros(szi), xy=zeros(szi), yx=zeros(szi), yy=zeros(szi),zx=zeros(szi),zy=zeros(szi)),
    #          j=(xx=zeros(szj), xy=zeros(szj), yx=zeros(szj), yy=zeros(szj),zx=zeros(szj),zy=zeros(szj)))
    # ε̇     = (i=(xx=zeros(szi), yy=zeros(szi), zz=zeros(szi), xy=zeros(szi),xz=zeros(szi),yz=zeros(szi)),
    #          j=(xx=zeros(szj), yy=zeros(szj), zz=zeros(szj), xy=zeros(szj),xz=zeros(szj),yz=zeros(szj))) 
    ε̇      = (xx=(i=zeros(szi),j=zeros(szj)) , yy=(i=zeros(szi),j=zeros(szj)),
              zz=(i=zeros(szi),j=zeros(szj)) , xy=(i=zeros(szi),j=zeros(szj)),
              xz=(i=zeros(szi),j=zeros(szj)) , yz=(i=zeros(szi),j=zeros(szj)))
    L      = (xx=(i=zeros(szi),j=zeros(szj)) , yy=(i=zeros(szi),j=zeros(szj)), zz=(i=zeros(szi), j=zeros(szj)),
              xy=(i=zeros(szi),j=zeros(szj)) , yx=(i=zeros(szi),j=zeros(szj)),
              zx=(i=zeros(szi),j=zeros(szj)) , zy=(i=zeros(szi),j=zeros(szj)))

    τ      = (xx=(i=zeros(szi),j=zeros(szj)),yy=(i=zeros(szi),j=zeros(szj)) , zz=(i=zeros(szi), j=zeros(szj)),
              xy=(i=zeros(szi),j=zeros(szj)),xz=(i=zeros(szi),j=zeros(szj)) , yz=(i=zeros(szi) ,j=zeros(szj)))
    
    τ0     = (xx=(i=zeros(szi),j=zeros(szj)),yy=(i=zeros(szi),j=zeros(szj)) , zz=(i=zeros(szi), j=zeros(szj)),
              xy=(i=zeros(szi),j=zeros(szj)),xz=(i=zeros(szi),j=zeros(szj)) , yz=(i=zeros(szi) ,j=zeros(szj)))
    # Storage on v and c meshes
    V     = ( x=(v=zeros(szv),c=zeros(szc)),
              y=(v=zeros(szv),c=zeros(szc)),
              z=(v=zeros(szv),c=zeros(szc))) 
    ρ     = (v=ones(szv)*ρ₀, c=ones(szc)*ρ₀)
    f_ext = (v=zeros(szv)  , c=zeros(szc))
    # BC
    Lbc        = 3.
    # BC on v and c mesh
    Δ0 = (x=0,y=0)
    Δi = (x=0,y=Δ.y/2)
    Δj = (x=Δ.x/2,y=0)
    bc_filt_V   = (v=Cerjean2D(X.x.v,X.y.v,Lbc,l,Δ0),c=Cerjean2D(X.x.c,X.y.c,Lbc,l,Δ))
    bc_filt_tau = (i=Cerjean2D(X.x.i,X.y.i,Lbc,l,Δi),j=Cerjean2D(X.x.j,X.y.j,Lbc,l,Δj))

    # Compute Ricker function with 2D spatial support
    f_ext = (v=zeros(szv)  , c=zeros(szc))
    xc2d   = X.x.c * ones(size( X.y.c))'
    yc2d   = ones(size( X.x.c)) * X.y.c'
    xv2d   = X.x.v * ones(size( X.y.v))'
    yv2d   = ones(size( X.x.v)) * X.y.v'

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
        for grid=1:2
            P0[grid] .= P[grid]]
        end

        for grid=1:2 ,comp=1:6 
            τ0[comp][grid] .= τ[comp][grid]
        end
    
        # 2D Ricker with spatial support
        @.. f_ext.c = ρ.c*Ricker.( xc2d, x₀, yc2d, y₀, t, t₀, 𝑓₀, σ₀)
        @.. f_ext.v = ρ.v*Ricker.( xv2d, x₀, yv2d, y₀, t, t₀, 𝑓₀, σ₀)
        
        # Inherited pressure (remove the instantaneous viscous contribution )
        for grid=1:2
             @.. P0[grid] = P0[grid] + χb(vol[grid]...)*∇V[grid]
        end

        # Inherited deviatoric stress (remove the instantaneous viscous contribution)
        for grid=1:2, comp=1:6 
            @.. τ0[comp][grid]= τ0[comp][grid] - χs(dev[grid]...)*ε̇[comp][grid]
        end
        

        # Velocity gradient components
        
        @.. L.xx.i[:,2:end-1] = (V.x.c[2:end,2:end-1] - V.x.c[1:end-1,2:end-1])/Δ.x
        @.. L.yx.i[:,2:end-1] = (V.y.c[2:end,2:end-1] - V.y.c[1:end-1,2:end-1])/Δ.x
        @.. L.zx.i[:,2:end-1] = (V.z.c[2:end,2:end-1] - V.z.c[1:end-1,2:end-1])/Δ.x

        @.. L.xx.j[2:end-1,:] = (V.x.v[2:end,:] - V.x.v[1:end-1,:])/Δ.x
        @.. L.yx.j[2:end-1,:] = (V.y.v[2:end,:] - V.y.v[1:end-1,:])/Δ.x
        @.. L.zx.j[2:end-1,:] = (V.z.v[2:end,:] - V.z.v[1:end-1,:])/Δ.x
  
        @.. L.yy.i[:,2:end-1] = (V.y.v[:,2:end] - V.y.v[:,1:end-1])/Δ.y
        @.. L.xy.i[:,2:end-1] = (V.x.v[:,2:end] - V.x.v[:,1:end-1])/Δ.y
        @.. L.zy.i[:,2:end-1] = (V.z.v[:,2:end] - V.z.v[:,1:end-1])/Δ.y

        @.. L.yy.j[2:end-1,:] = (V.y.c[2:end-1,2:end] - V.y.c[2:end-1,1:end-1])/Δ.y
        @.. L.xy.j[2:end-1,:] = (V.x.c[2:end-1,2:end] - V.x.c[2:end-1,1:end-1])/Δ.y
        @.. L.zy.j[2:end-1,:] = (V.z.c[2:end-1,2:end] - V.z.c[2:end-1,1:end-1])/Δ.y
        

        
        for grid=1:2
        L.xy[grid][:,topbc[grid]] .= (-L.yx[grid][:,topbc[grid]] .* ηs(dev[grid]...)[:,topbc[grid]] .- 
                                       θs(dev[grid]...)[:,topbc[grid]].*τ0.xy[grid][:,topbc[grid]]) ./ 
                                       ηs(dev[grid]...)[:,topbc[grid]]
        L.zy[grid][:,topbc[grid]] .= -θs(dev[grid]...)[:,topbc[grid]].*  τ0.yz[grid][:,topbc[grid]] ./
                                       ηs(dev[grid]...)[:,topbc[grid]]
        L.yy[grid][:,topbc[grid]] .= (- 3.0 * L.xx[grid][:,topbc[grid]] .* ηb(vol[grid]...)[:,topbc[grid]] + 
                                        2.0 * L.xx[grid][:,topbc[grid]] .* ηs(dev[grid]...)[:,topbc[grid]] +
                                        3.0 * P0[grid][:,topbc[grid]]   .* θb(vol[grid]...)[:,topbc[grid]] - 
                                        3.0 * θs(dev[grid]...)[:,topbc[grid]]  .* τ0.yy[grid][:,topbc[grid]]) ./
                                       (3.0 * ηb(vol[grid]...)[:,topbc[grid]]  + 4.0 * ηs(dev[grid]...)[:,topbc[grid]])
        end
    
        # Divergence
        for grid=1:2
            @.. ∇V[grid]   = L.xx[grid] + L.yy[grid]
        end

        # Deviatoric strain rate 
        for grid=1:2
            @.. ε̇.xx[grid] = L.xx[grid] - 1//3*∇V[grid]
            @.. ε̇.yy[grid] = L.yy[grid] - 1//3*∇V[grid]
            @.. ε̇.zz[grid] = - 1//3*∇V[grid]
            @.. ε̇.xy[grid] = 1//2*(L.xy[grid] + L.yx[grid])
            # in 2D Lxz and Lyz are zero 
            @.. ε̇.xz[grid] = 1//2*(L.zx[grid])
            @.. ε̇.yz[grid] = 1//2*(L.zy[grid])
        end

        # Stress update

        for grid=1:2, comp=1:6
             @.. τ[comp][grid] = ηs(dev[grid]...)*(ε̇[comp][grid]) + θs(dev[grid]...)*τ0[comp][grid]
        end        
        # Pressure update 
        for grid=1:2
            @.. P[grid]  = θb(vol[grid]...)*P0[grid] - ηb(vol[grid]...)*∇V[grid]
        end 
        #@.. P.j    = θb(vol.j...)*P0.j - ηb(vol.j...)*∇V.j 

        # these stresses are on the free surface) 
        # τ.j.xy[:,end] .= 0.
        # τ.j.yz[:,end] .= 0.
        # τ.j.yy[:,end] .= P.j[:,end]

        # Linear momentum balance
        @.. V.x.v[2:end-1,2:end-1] = (V.x.v[2:end-1,2:end-1] 
                                    + Δt/ρ.v[2:end-1,2:end-1]
                                    *((τ.xx.j[3:end-1,2:end-1]-τ.xx.j[2:end-2,2:end-1])/Δ.x
                                    + (τ.xy.i[2:end-1,3:end-1]-τ.xy.i[2:end-1,2:end-2])/Δ.y 
                                    - (P.j[3:end-1,2:end-1]-P.j[2:end-2,2:end-1])/Δ.x 
                                    - facS.v.x*f_ext.v[2:end-1,2:end-1]))

        @.. V.x.c[2:end-1,2:end-1] = (V.x.c[2:end-1,2:end-1] 
                                    + Δt/ρ.c[2:end-1,2:end-1]
                                    *((τ.xx.i[2:end,2:end-1]-τ.xx.i[1:end-1,2:end-1])/Δ.x
                                    + (τ.xy.j[2:end-1,2:end]-τ.xy.j[2:end-1,1:end-1])/Δ.y
                                    - (P.i[2:end,2:end-1]-P.i[1:end-1,2:end-1])/Δ.x 
                                    - facS.c.x*f_ext.c[2:end-1,2:end-1]))                            

        @.. V.y.v[2:end-1,2:end-1] = (V.y.v[2:end-1,2:end-1] 
                                    + Δt/ρ.v[2:end-1,2:end-1]
                                    *((τ.xy.j[3:end-1,2:end-1]-τ.xy.j[2:end-2,2:end-1])/Δ.x
                                    + (τ.yy.i[2:end-1,3:end-1]-τ.yy.i[2:end-1,2:end-2])/Δ.y 
                                    - (P.i[2:end-1,3:end-1]-P.i[2:end-1,2:end-2])/Δ.y 
                                    - facS.v.y*f_ext.v[2:end-1,2:end-1]))
        
        @.. V.y.c[2:end-1,2:end-1] = (V.y.c[2:end-1,2:end-1] 
                                    + Δt/ρ.c[2:end-1,2:end-1]
                                    *((τ.xy.i[2:end,2:end-1]-τ.xy.i[1:end-1,2:end-1])/Δ.x
                                    + (τ.yy.j[2:end-1,2:end]-τ.yy.j[2:end-1,1:end-1])/Δ.y 
                                    - (P.j[2:end-1,2:end]-P.j[2:end-1,1:end-1])/Δ.y 
                                    - facS.c.y*f_ext.c[2:end-1,2:end-1]))   

        # the two terms in dPdz and dtauzzdz  cancel in linear elastic case ... but i am not sure with other rheologies so I have left them 
        @.. V.z.v[2:end-1,2:end-1] = (V.z.v[2:end-1,2:end-1] 
                                    + Δt/ρ.v[2:end-1,2:end-1]
                                    *((τ.xz.j[3:end-1,2:end-1]-τ.xz.j[2:end-2,2:end-1])/Δ.x
                                    + (τ.yz.i[2:end-1,3:end-1]-τ.yz.i[2:end-1,2:end-2])/Δ.y 
                                    - facS.v.z* f_ext.v[2:end-1,2:end-1]))
        
        @.. V.z.c[2:end-1,2:end-1] = (V.z.c[2:end-1,2:end-1] 
                                    + Δt/ρ.c[2:end-1,2:end-1]
                                    *((τ.xz.i[2:end,2:end-1]-τ.xz.i[1:end-1,2:end-1])/Δ.x
                                    + (τ.yz.j[2:end-1,2:end]-τ.yz.j[2:end-1,1:end-1])/Δ.y 
                                    - facS.c.z*f_ext.c[2:end-1,2:end-1]))  
                                    
        # free surface update
        @.. V.x.v[2:end-1,end] = (V.x.v[2:end-1,end] 
                                 + Δt/ρ.v[2:end-1,end]
                                 *((τ.xx.j[3:end-1,end]-τ.xx.j[2:end-2,end])/Δ.x
                                 + (τ.xy.i[2:end-1,end]-τ.xy.i[2:end-1,end-1])/Δ.y 
                                 - (P.j[3:end-1,end]-P.j[2:end-2,end])/Δ.x 
                                 - facS.v.x*f_ext.v[2:end-1,end]))

        @.. V.y.v[2:end-1,end] = (V.y.v[2:end-1,end] 
                                 + Δt/ρ.v[2:end-1,end]
                                 *((τ.xy.j[3:end-1,end]-τ.xy.j[2:end-2,end])/Δ.x
                                 + (τ.yy.i[2:end-1,end]-τ.yy.i[2:end-2,end-1])/Δ.y 
                                 - (P.i[2:end-1,end]   -P.i[2:end-1,end-1])/Δ.y 
                                 - facS.v.y*f_ext.v[2:end-1,end]))

        @.. V.z.v[2:end-1,end] = (V.z.v[2:end-1,end] 
                                    + Δt/ρ.v[2:end-1,end]
                                    *((τ.xz.j[3:end-1,end]-τ.xz.j[2:end-2,end])/Δ.x
                                    + (τ.yz.i[2:end-1,end]-τ.yz.i[2:end-1,end-1])/Δ.y 
                                    - facS.v.z* f_ext.v[2:end-1,end]))


    
        # Absorbing boundary Cerjean et al. (1985)
        for grid=1:2, comp=1:length(Δ) 
            @..  V[comp][grid] = V[comp][grid]  * bc_filt_V[grid] 
        end 
 
        for grid=1:2
            @..  P[grid]    = P[grid]    *  bc_filt_tau[grid] 
        end
        
        for grid=1:2, comp=1:6 
            @..  τ[comp][grid] *= bc_filt_V[grid] 
        end 


        # Visualisation
        if mod(it, Nout)==0 && visu==true

            resol=200 
            f = Figure(resolution = (l.x/l.y*resol*2, resol*2), fontsize=15)
            for grid=1:2, comp=1:length(Δ) 
                ax = Axis(f[grid, comp], aspect=l.x/l.y, title = L" v%$(comp) on grid %$(grid) ", xlabel = L"$x$ [m]", ylabel = L"$y$ [m]")
                hm = GLMakie.heatmap!(ax, X.x.v, X.y.v, V.x.v, colormap = wave_colors,colorrange=(-3.e-5,3.e-5))
            end

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
         velocity_matrix.x[:,ksismo] .= V.x.v[:,end-1]
         velocity_matrix.y[:,ksismo] .= V.y.v[:,end-1]
         velocity_matrix.z[:,ksismo] .= V.z.v[:,end-1]
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

function Cerjean2D(x,y,Lbc,l,Δ)
    return ((1.0 .- exp.(-(x*ones(size(y))'.- (0*l.x+Δ.x/2)).^2/Lbc.^2))
         .*(1.0 .- exp.(-(x*ones(size(y))' .- (l.x-Δ.x/2)).^2/Lbc.^2))
         .*(1.0 .- exp.(-(ones(size(x))*y' .- (0*l.y+Δ.y/2)).^2/Lbc.^2))
         )
         #.*(1.0 .- exp.(-(ones(size(X.x))*X.y' .-  l.y).^2/Lbc.^2))
end

function Print2Disk( f, path, field, istep; res=4)
    path1 = path*"/_$field/"
    mkpath(path1)
    save(path1*"$field"*@sprintf("%05d", istep)*".png", f, px_per_unit = res) 
end

MainSource()
using SeismicQ, FastBroadcast, GLMakie, Printf, Colors, ColorSchemes, MathTeXEngine, UnPack, Makie.GeometryBasics
Makie.update_theme!(fonts = (regular = texfont(), bold = texfont(:bold), italic = texfont(:italic)))

@views h(x,A,σ,b,x0)    = A*exp(-(x-x0)^2/σ^2) + b
@views dhdx(x,A,σ,b,x0) = -2*x/σ^2*A*exp(-(x-x0).^2/σ^2)


function Mesh_y( X, A, x0, σ, b, m, ymin0, ymax0, σy )
    y0    = ymax0
    ymin1 = (sinh.( σy.*(ymin0.-y0) ))
    ymax1 = (sinh.( σy.*(ymax0.-y0) ))
    sy    = (ymax0-ymin0)/(ymax1-ymin1)
    y     = (sinh.( σy.*(X[2].-y0) )) .* sy  .+ y0
    # y = X[2]
    z0    = -(A*exp(-(X[1]-x0)^2/σ^2) + b) # topography height
    y     = (y/ymin0)*((z0+m))-z0        # shift grid vertically
    return y
end
function Mesh_x( X, A, x0, σ, b, m, xmin0, xmax0, σx )
    xmin1 = (sinh.( σx.*(xmin0.-x0) ))
    xmax1 = (sinh.( σx.*(xmax0.-x0) ))
    sx    = (xmax0-xmin0)/(xmax1-xmin1)
    x     = (sinh.( σx.*(X[1].-x0) )) .* sx  .+ x0        
    # x   = X[1]
    return x
end

function PatchPlotMakie(vertx, verty, field; cmap = :turbo, write_fig=false )
    f    = GLMakie.Figure(resolution = (1200, 1000))
    xmin = minimum(vertx)
    xmax = maximum(vertx)
    ymin = minimum(verty)
    ymax = maximum(verty)
    ar = (xmax - xmin) / (ymax - ymin)
    GLMakie.Axis(f[1,1]) #, aspect = ar
    min_v = minimum( field ); max_v = maximum( field )
    limits = min_v ≈ max_v ? (min_v, min_v + 1) : (min_v, max_v)
    p = [Polygon( Point2f0[ (vertx[i,j], verty[i,j]) for j=1:4] ) for i in 1:length(field)]
    GLMakie.poly!(p, color = field, colormap = cmap, strokewidth = 1, strokecolor = :black, markerstrokewidth = 0, markerstrokecolor = (0, 0, 0, 0), aspect=:image, colorrange=limits)
    GLMakie.Colorbar(f[1, 2], colormap = cmap, limits=limits, flipaxis = true, size = 25 )
    display(f)
    # if write_fig==true 
    #     FileIO.save( string(@__DIR__, "/plot.png"), f)
    # end
    return nothing
end

function InverseJacobian!(∂ξ,∂η,∂x,∂y)
    M = zeros(2,2)
    @time for i in eachindex(∂ξ.∂x)
        M[1,1]   = ∂x.∂ξ[i]
        M[1,2]   = ∂x.∂η[i]
        M[2,1]   = ∂y.∂ξ[i]
        M[2,2]   = ∂y.∂η[i]
        invJ     = inv(M)
        ∂ξ.∂x[i] = invJ[1,1]
        ∂ξ.∂y[i] = invJ[1,2]
        ∂η.∂x[i] = invJ[2,1]
        ∂η.∂y[i] = invJ[2,2]
    end
    @printf("min(∂ξ∂x) = %1.6f --- max(∂ξ∂x) = %1.6f\n", minimum(∂ξ.∂x), maximum(∂ξ.∂x))
    @printf("min(∂ξ∂y) = %1.6f --- max(∂ξ∂y) = %1.6f\n", minimum(∂ξ.∂y), maximum(∂ξ.∂y))
    @printf("min(∂η∂x) = %1.6f --- max(∂η∂x) = %1.6f\n", minimum(∂η.∂x), maximum(∂η.∂x))
    @printf("min(∂η∂y) = %1.6f --- max(∂η∂y) = %1.6f\n", minimum(∂η.∂y), maximum(∂η.∂y))
    return nothing
end

function ComputeForwardTransformation_ini!( ∂x, ∂y, x_ini, y_ini, X_msh, Amp, x0, σ, m, x, y, σx, σy, ϵ)
    
    xmin, xmax = x.min, x.max
    ymin, ymax = y.min, y.max

    @time for i in eachindex(y_ini)          
    
        # compute dxdksi
        X_msh[1] = x_ini[i]-ϵ
        X_msh[2] = y_ini[i] 
        xm       = Mesh_x( X_msh,  Amp, x0, σ, xmax, m, xmin, xmax, σx )
        # --------
        X_msh[1] = x_ini[i]+ϵ
        X_msh[2] = y_ini[i]
        xp       = Mesh_x( X_msh,  Amp, x0, σ, xmax, m, xmin, xmax, σx )
        # --------
        ∂x.∂ξ[i] = (xp - xm) / (2ϵ)
    
        # compute dydeta
        X_msh[1] = x_ini[i]
        X_msh[2] = y_ini[i]-ϵ
        xm     = Mesh_x( X_msh,  Amp, x0, σ, ymax, m, ymin, ymax, σy )
        # --------
        X_msh[1] = x_ini[i]
        X_msh[2] = y_ini[i]+ϵ
        xp       = Mesh_x( X_msh,  Amp, x0, σ, ymax, m, ymin, ymax, σy )
        # --------
        ∂x.∂η[i] = (xp - xm) / (2ϵ)
    
        # compute dydksi
        X_msh[1] = x_ini[i]-ϵ
        X_msh[2] = y_ini[i] 
        ym       = Mesh_y( X_msh,  Amp, x0, σ, ymax, m, ymin, ymax, σy )
        # --------
        X_msh[1] = x_ini[i]+ϵ
        X_msh[2] = y_ini[i]
        yp       = Mesh_y( X_msh,  Amp, x0, σ, ymax, m, ymin, ymax, σy )
        # --------
        ∂y.∂ξ[i] = (yp - ym) / (2ϵ)
    
        # compute dydeta
        X_msh[1] = x_ini[i]
        X_msh[2] = y_ini[i]-ϵ
        ym     = Mesh_y( X_msh,  Amp, x0, σ, ymax, m, ymin, ymax, σy )
        # --------
        X_msh[1] = x_ini[i]
        X_msh[2] = y_ini[i]+ϵ
        yp     = Mesh_y( X_msh,  Amp, x0, σ, ymax, m, ymin, ymax, σy )
        # --------
        ∂y.∂η[i] = (yp - ym) / (2ϵ)
    end
    # #################
    # # ForwardDiff
    # g = zeros(2)
    # Y = zeros(1)
    # dydksi_FD = zeros(size(dydeta))
    # dydeta_FD = zeros(size(dydeta))
    # dxdksi_FD = zeros(size(dydeta))
    # dxdeta_FD = zeros(size(dydeta))
    # @time for i in eachindex(dydeta_FD)
    #     X_msh[1] = x_ini[i]
    #     X_msh[2] = y_ini[i]
    #     Mesh_y_closed = (X_msh) -> Mesh_y( X_msh, Amp, x0, σ, b, m, ymin )
    #     ForwardDiff.gradient!( g, Mesh_y_closed, X_msh )
    #     dydksi_FD[i] = g[1]
    #     dydeta_FD[i] = g[2]
    #     Meshx_surf_closed = (X_msh) -> Mesh_x( X_msh, Amp, x0, σ, b, m, ymin )
    #     ForwardDiff.gradient!( g, Meshx_surf_closed, X_msh )
    #     dxdksi_FD[i] = g[1]
    #     dxdeta_FD[i] = g[2]
    # end
    
    # dxdksi_num = diff(xv4,dims=1)/(Δx/2)
    # dxdeta_num = diff(xv4,dims=2)/(Δy/2)
    # dydksi_num = diff(yv4,dims=1)/(Δx/2)
    # dydeta_num = diff(yv4,dims=2)/(Δy/2)
    
    # @printf("min(dxdksi    ) = %1.6f --- max(dxdksi    ) = %1.6f\n", minimum(dxdksi   ), maximum(dxdksi   ))
    # @printf("min(dxdksi_FD ) = %1.6f --- max(dxdksi_FD ) = %1.6f\n", minimum(dxdksi_FD), maximum(dxdksi_FD))
    # @printf("min(dxdksi_num) = %1.6f --- max(dxdksi_num) = %1.6f\n", minimum(dxdksi_num), maximum(dxdksi_num))
    
    # @printf("min(dxdeta    ) = %1.6f --- max(dxdeta   ) = %1.6f\n", minimum(dxdeta   ), maximum(dxdeta   ))
    # @printf("min(dxdeta_FD ) = %1.6f --- max(dxdeta_FD) = %1.6f\n", minimum(dxdeta_FD), maximum(dxdeta_FD))
    # @printf("min(dxdeta_num) = %1.6f --- max(dxdeta_num) = %1.6f\n", minimum(dxdeta_num), maximum(dxdeta_num))
    
    # @printf("min(dydksi    ) = %1.6f --- max(dydksi    ) = %1.6f\n", minimum(dydksi   ), maximum(dydksi   ))
    # @printf("min(dydksi_FD ) = %1.6f --- max(dydksi_FD ) = %1.6f\n", minimum(dydksi_FD), maximum(dydksi_FD))
    # @printf("min(dydksi_num) = %1.6f --- max(dydksi_num) = %1.6f\n", minimum(dydksi_num), maximum(dydksi_num))
    
    # @printf("min(dydeta    ) = %1.6f --- max(dydeta    ) = %1.6f\n", minimum(dydeta   ), maximum(dydeta   ))
    # @printf("min(dydeta_FD ) = %1.6f --- max(dydeta_FD ) = %1.6f\n", minimum(dydeta_FD), maximum(dydeta_FD))
    # @printf("min(dydeta_num) = %1.6f --- max(dydeta_num) = %1.6f\n", minimum(dydeta_num), maximum(dydeta_num))
    return nothing
end

function MainSource()
    adapt_mesh = true
    visu     = true
    printfig = false  # print figures to disk
    path     = "./runs/"
    juliadivcmap    = zeros(RGB{Float64}, 5)
    juliadivcmap[1] = RGBA{Float64}(0/255,150/255,0/255, 1.)  
    juliadivcmap[2] = RGBA{Float64}(0/255,0/255,200/255, 1.)  
    juliadivcmap[3] = RGBA{Float64}(255/255,255/255,255/255, 1.) 
    juliadivcmap[4] = RGBA{Float64}(150/255,0/255,150/255, 1.) 
    juliadivcmap[5] = RGBA{Float64}(200/255,0/255,0/255, 1.)
    wave_colors     = cgrad(juliadivcmap, length(juliadivcmap), categorical=false, rev=false)
    
    # Spatial extent
    l  = (x = 25, y = 25)
    x  = (min=-l.x/2, max=l.x/2)
    y  = (min=-l.y/1, max=0.)

    # Discretization
    Nc  = (x = 50, y = 50) 
    Δ   = (ξ = l.x/Nc.x, η = l.y/Nc.y, ζ = 1.0)
    X   = (v = (x = LinRange(x.min,       x.max,       Nc.x+1) , y = LinRange(y.min,       y.max,Nc.y+1)),
           c = (x = LinRange(x.min-Δ.ξ/2, x.max+Δ.ξ/2, Nc.x+2) , y = LinRange(y.min-Δ.η/2, y.max+Δ.η/2,Nc.y+2)),
           i = (x = LinRange(x.min,       x.max,       Nc.x+1) , y = LinRange(y.min-Δ.η/2, y.max+Δ.η/2,Nc.y+2)),
           j = (x = LinRange(x.min-Δ.ξ/2, x.max+Δ.ξ/2, Nc.x+2) , y = LinRange(y.min,       y.max,Nc.y+1))) 
        
    # Source parameters
    𝑓₀   = 100   # Central frequency of the source [Hz]
    t₀   = 1.2/𝑓₀
    σ₀   = l.x/100
    x₀   = (x.min + x.max)/2
    y₀   = (y.min + y.max)/2
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
    VolRheo = :KelvinVE  #:Elastic or :KelvinVE 
    
    DevRheo = :Elastic #or :MaxwellVE
    #VolRheo = :Elastic #or :KelvinVE 

    # Time domain
    c_eff = sqrt((K₀*(1+Fb_b)+4/3*G₀)/ρ₀) 
    Δt    = min(1e10, 0.1*Δ.ξ/c_eff, 0.1*Δ.η/c_eff ) # Courant criteria from wavespeed
    Nt    = 1000
    Nout  = 1000
    t     = -t₀
   
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
    P0    = (i=zeros(szi),     j=zeros(szj))
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
    f_ext  = (v=zeros(szv)  , c=zeros(szc))
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
    # Initialisation
    xxv, yyv    = LinRange(x.min-Δ.ξ/2, x.max+Δ.ξ/2, 2Nc.x+3), LinRange(y.min-Δ.η/2, y.max+Δ.η/2, 2Nc.y+3)
    (xv4,yv4) = ([x for x=xxv,y=yyv], [y for x=xxv,y=yyv])
    ∂ξ∂x =  ones(2Nc.x+3, 2Nc.y+3)
    ∂ξ∂y = zeros(2Nc.x+3, 2Nc.y+3)
    ∂η∂x = zeros(2Nc.x+3, 2Nc.y+3)
    ∂η∂y =  ones(2Nc.x+3, 2Nc.y+3)
    hx   = zeros(2Nc.x+3, 2Nc.y+3)
    if adapt_mesh
        x0     = (x.min + x.max)/2
        m      = y.min
        Amp    = 2.0
        σ      = 0.9
        σx     = 0.1
        σy     = 0.1
        ϵ      = 1e-7
        # copy initial y
        x_ini  = copy(xv4)
        y_ini  = copy(yv4)
        X_msh  = zeros(2)
        # Compute slope
        hx     = -dhdx.(x_ini, Amp, σ, y.max, x0)
        # Deform mesh
        for i in eachindex(x_ini)          
            X_msh[1] = x_ini[i]
            X_msh[2] = y_ini[i]     
            xv4[i]   =  Mesh_x( X_msh,  Amp, x0, σ, y.max, m, x.min, x.max, σx )
            yv4[i]   =  Mesh_y( X_msh,  Amp, x0, σ, y.max, m, y.min, y.max, σy )
        end
        # Compute forward transformation
        # params = (Amp=Amp, x0=x0, σ=σ, m=m, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, σx=σx, σy=σy, ϵ=ϵ)
        ∂x     = (∂ξ=zeros(size(yv4)), ∂η = zeros(size(yv4)) )
        ∂y     = (∂ξ=zeros(size(yv4)), ∂η = zeros(size(yv4)) )
        ComputeForwardTransformation_ini!( ∂x, ∂y, x_ini, y_ini, X_msh, Amp, x0, σ, m, x, y, σx, σy, ϵ)
        # Solve for inverse transformation
        ∂ξ = (∂x=∂ξ∂x, ∂y=∂ξ∂y); ∂η = (∂x=∂η∂x, ∂y=∂η∂y)
        InverseJacobian!(∂ξ,∂η,∂x,∂y)
        ∂ξ∂x .= ∂ξ.∂x; ∂ξ∂y .= ∂ξ.∂y
        ∂η∂x .= ∂η.∂x; ∂η∂y .= ∂η.∂y
    end
    ∂ξ∂x1 = (
        i = ∂ξ∂x[3:2:end-2,2:2:end-1],
        j = ∂ξ∂x[3:2:end-2,2:2:end-1],
        c = ∂ξ∂x[1:2:end-0,1:2:end-0],
        v = ∂ξ∂x[2:2:end-1,2:2:end-1],
    )  
    ∂η∂x1 = (
        i = ∂η∂x[3:2:end-2,2:2:end-1],
        j = ∂η∂x[3:2:end-2,2:2:end-1],
        c = ∂η∂x[1:2:end-0,1:2:end-0],
        v = ∂η∂x[2:2:end-1,2:2:end-1],
    ) 
    ∂ξ∂y1 = (
        i = ∂ξ∂y[3:2:end-2,2:2:end-1],
        j = ∂ξ∂y[3:2:end-2,2:2:end-1],
        c = ∂ξ∂y[1:2:end-0,1:2:end-0],
        v = ∂ξ∂y[2:2:end-1,2:2:end-1],
    )  
    ∂η∂y1 = (
        i = ∂η∂y[3:2:end-2,2:2:end-1],
        j = ∂η∂y[3:2:end-2,2:2:end-1],
        c = ∂η∂y[1:2:end-0,1:2:end-0],
        v = ∂η∂y[2:2:end-1,2:2:end-1],
    )
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
        @.. L.i.xx[:,2:end-1] = ∂ξ∂x1.i * (V.c.x[2:end,2:end-1] - V.c.x[1:end-1,2:end-1])/Δ.ξ + ∂η∂x1.i * (V.v.x[ :     ,2:end] - V.v.x[ :     ,1:end-1])/Δ.η
        @.. L.j.xx[2:end-1,:] = ∂ξ∂x1.j * (V.v.x[2:end, :     ] - V.v.x[1:end-1, :     ])/Δ.ξ + ∂η∂x1.j * (V.c.x[2:end-1,2:end] - V.c.x[2:end-1,1:end-1])/Δ.η

        @.. L.i.yx[:,2:end-1] = ∂ξ∂x1.i * (V.c.y[2:end,2:end-1] - V.c.y[1:end-1,2:end-1])/Δ.ξ + ∂η∂x1.i * (V.v.y[ :     ,2:end] - V.v.y[ :     ,1:end-1])/Δ.η
        @.. L.j.yx[2:end-1,:] = ∂ξ∂x1.j * (V.v.y[2:end, :     ] - V.v.y[1:end-1, :     ])/Δ.ξ + ∂η∂x1.j * (V.c.y[2:end-1,2:end] - V.c.y[2:end-1,1:end-1])/Δ.η

        @.. L.i.yy[:,2:end-1] = ∂ξ∂y1.i * (V.c.y[2:end,2:end-1] - V.c.y[1:end-1,2:end-1])/Δ.ξ + ∂η∂y1.i * (V.v.y[ :     ,2:end] - V.v.y[ :     ,1:end-1])/Δ.η
        @.. L.j.yy[2:end-1,:] = ∂ξ∂y1.j * (V.v.y[2:end, :     ] - V.v.y[1:end-1, :     ])/Δ.ξ + ∂η∂y1.j * (V.c.y[2:end-1,2:end] - V.c.y[2:end-1,1:end-1])/Δ.η

        @.. L.i.xy[:,2:end-1] = ∂ξ∂y1.i * (V.c.x[2:end,2:end-1] - V.c.x[1:end-1,2:end-1])/Δ.ξ + ∂η∂y1.i * (V.v.x[ :,     2:end] - V.v.x[ :     ,1:end-1])/Δ.η
        @.. L.j.xy[2:end-1,:] = ∂ξ∂y1.j * (V.v.x[2:end, :     ] - V.v.x[1:end-1, :     ])/Δ.ξ + ∂η∂y1.j * (V.c.x[2:end-1,2:end] - V.c.x[2:end-1,1:end-1])/Δ.η

        @.. L.i.zy[:,2:end-1] = (V.v.z[ :     ,2:end] - V.v.z[ :     ,1:end-1])/Δ.η
        @.. L.j.zy[2:end-1,:] = (V.c.z[2:end-1,2:end] - V.c.z[2:end-1,1:end-1])/Δ.η

        @.. L.i.zx[:,2:end-1] = (V.c.z[2:end,2:end-1] - V.c.z[1:end-1,2:end-1])/Δ.ξ
        @.. L.j.zx[2:end-1,:] = (V.v.z[2:end, :     ] - V.v.z[1:end-1, :     ])/Δ.ξ
        
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

        # @show size(τ.j.xx[3:end-1,2:end-1])
        # @show size(V.v.x[2:end-1,2:end-1])
        # @show size(∂ξ∂x1.v[2:end-1,2:end-1])

        # @show size(τ.i.xx[2:end,2:end-1])
        # @show size(V.c.x[2:end-1,2:end-1])
        # @show size(∂ξ∂x1.c[2:end-1,2:end-1])

        # Linear momentum balance
        @. V.v.x[2:end-1,2:end-1] = (V.v.x[2:end-1,2:end-1] 
                                    + Δt/ρ.v[2:end-1,2:end-1]
                                    *( ∂ξ∂x1.v[2:end-1,2:end-1] * (τ.j.xx[3:end-1,2:end-1]-τ.j.xx[2:end-2,2:end-1])/Δ.ξ + ∂η∂x1.v[2:end-1,2:end-1] * (τ.i.xx[2:end-1,3:end-1]-τ.i.xx[2:end-1,2:end-2])/Δ.η
                                    +  ∂ξ∂y1.v[2:end-1,2:end-1] * (τ.j.xy[3:end-1,2:end-1]-τ.j.xy[2:end-2,2:end-1])/Δ.ξ + ∂η∂y1.v[2:end-1,2:end-1] * (τ.i.xy[2:end-1,3:end-1]-τ.i.xy[2:end-1,2:end-2])/Δ.η 
                                    -  ∂ξ∂x1.v[2:end-1,2:end-1] * (   P.j[3:end-1,2:end-1]-   P.j[2:end-2,2:end-1])/Δ.ξ - ∂η∂x1.v[2:end-1,2:end-1] * (   P.i[2:end-1,3:end-1]-   P.i[2:end-1,2:end-2])/Δ.η  
                                    -  facS.v.x*f_ext.v[2:end-1,2:end-1]))
        @. V.c.x[2:end-1,2:end-1] = (V.c.x[2:end-1,2:end-1] 
                                    + Δt/ρ.c[2:end-1,2:end-1]
                                    *( ∂ξ∂x1.c[2:end-1,2:end-1] * (τ.i.xx[2:end,2:end-1]-τ.i.xx[1:end-1,2:end-1])/Δ.ξ + ∂η∂x1.c[2:end-1,2:end-1] * (τ.j.xx[2:end-1,2:end]-τ.j.xx[2:end-1,1:end-1])/Δ.η 
                                    +  ∂ξ∂y1.c[2:end-1,2:end-1] * (τ.i.xy[2:end,2:end-1]-τ.i.xy[1:end-1,2:end-1])/Δ.ξ + ∂η∂y1.c[2:end-1,2:end-1] * (τ.j.xy[2:end-1,2:end]-τ.j.xy[2:end-1,1:end-1])/Δ.η
                                    -  ∂ξ∂x1.c[2:end-1,2:end-1] * (   P.i[2:end,2:end-1]-   P.i[1:end-1,2:end-1])/Δ.ξ + ∂η∂x1.c[2:end-1,2:end-1] * (   P.j[2:end-1,2:end]-   P.j[2:end-1,1:end-1])/Δ.η 
                                    -  facS.c.x*f_ext.c[2:end-1,2:end-1]))                            

        @. V.v.y[2:end-1,2:end-1] = (V.v.y[2:end-1,2:end-1] 
                                    + Δt/ρ.v[2:end-1,2:end-1]
                                    *( ∂ξ∂x1.v[2:end-1,2:end-1] * (τ.j.xy[3:end-1,2:end-1]-τ.j.xy[2:end-2,2:end-1])/Δ.ξ + ∂η∂x1.v[2:end-1,2:end-1] * (τ.i.xy[2:end-1,3:end-1]-τ.i.xy[2:end-1,2:end-2])/Δ.η 
                                    +  ∂ξ∂y1.v[2:end-1,2:end-1] * (τ.j.yy[3:end-1,2:end-1]-τ.j.yy[2:end-2,2:end-1])/Δ.ξ + ∂η∂y1.v[2:end-1,2:end-1] * (τ.i.yy[2:end-1,3:end-1]-τ.i.yy[2:end-1,2:end-2])/Δ.η 
                                    -  ∂ξ∂y1.v[2:end-1,2:end-1] * (   P.j[3:end-1,2:end-1]-   P.j[2:end-2,2:end-1])/Δ.ξ - ∂η∂y1.v[2:end-1,2:end-1] * (   P.i[2:end-1,3:end-1]-   P.i[2:end-1,2:end-2])/Δ.η 
                                    -  facS.v.y*f_ext.v[2:end-1,2:end-1]))
        
        @.. V.c.y[2:end-1,2:end-1] = (V.c.y[2:end-1,2:end-1] 
                                    + Δt/ρ.c[2:end-1,2:end-1]
                                    *( ∂ξ∂x1.c[2:end-1,2:end-1] * (τ.i.xy[2:end,2:end-1]-τ.i.xy[1:end-1,2:end-1])/Δ.ξ + ∂η∂x1.c[2:end-1,2:end-1] * (τ.j.xy[2:end-1,2:end]-τ.j.xy[2:end-1,1:end-1])/Δ.η
                                    +  ∂ξ∂y1.c[2:end-1,2:end-1] * (τ.i.yy[2:end,2:end-1]-τ.i.yy[1:end-1,2:end-1])/Δ.ξ + ∂η∂y1.c[2:end-1,2:end-1] * (τ.j.yy[2:end-1,2:end]-τ.j.yy[2:end-1,1:end-1])/Δ.η 
                                    -  ∂ξ∂y1.c[2:end-1,2:end-1] * (   P.i[2:end,2:end-1]-   P.i[1:end-1,2:end-1])/Δ.ξ - ∂η∂y1.c[2:end-1,2:end-1] * (P.j[2:end-1,2:end]-P.j[2:end-1,1:end-1])/Δ.η 
                                    -  facS.c.y*f_ext.c[2:end-1,2:end-1]))   

        # the two terms in dPdz and dtauzzdz  cancel in linear elastic case ... but i am not sure with other rheologies so I have left them 
        @.. V.v.z[2:end-1,2:end-1] = (V.v.z[2:end-1,2:end-1] 
                                    + Δt/ρ.v[2:end-1,2:end-1]
                                    *( (τ.j.xz[3:end-1,2:end-1]-τ.j.xz[2:end-2,2:end-1])/Δ.ξ
                                    +  (τ.i.yz[2:end-1,3:end-1]-τ.i.yz[2:end-1,2:end-2])/Δ.η 
                                    -  facS.v.z* f_ext.v[2:end-1,2:end-1]))
        
        @.. V.c.z[2:end-1,2:end-1] = (V.c.z[2:end-1,2:end-1] 
                                    + Δt/ρ.c[2:end-1,2:end-1]
                                    *( (τ.i.xz[2:end,2:end-1]-τ.i.xz[1:end-1,2:end-1])/Δ.ξ
                                    +  (τ.j.yz[2:end-1,2:end]-τ.j.yz[2:end-1,1:end-1])/Δ.η 
                                    -  facS.c.z*f_ext.c[2:end-1,2:end-1]))   
    
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
            xv2_1, yv2_1 = xv4[2:2:end-1,2:2:end-1  ], yv4[2:2:end-1,2:2:end-1  ]
            xv2_2, yv2_2 = xv4[1:2:end-0,1:2:end-0  ], yv4[1:2:end-0,1:2:end-0  ]
            xc2_1, yc2_1 = xv4[3:2:end-2,2:2:end-1  ], yv4[3:2:end-2,2:2:end-1  ]
            xc2_2, yc2_2 = xv4[2:2:end-1,3:2:end-2+2], yv4[2:2:end-1,3:2:end-2+2]
            vertx = [  xv2_1[1:end-1,1:end-1][:]  xv2_1[2:end-0,1:end-1][:]  xv2_1[2:end-0,2:end-0][:]  xv2_1[1:end-1,2:end-0][:] ] 
            verty = [  yv2_1[1:end-1,1:end-1][:]  yv2_1[2:end-0,1:end-1][:]  yv2_1[2:end-0,2:end-0][:]  yv2_1[1:end-1,2:end-0][:] ] 
            PatchPlotMakie(vertx, verty, V.c.x[2:end-1,2:end-1][:]; cmap=wave_colors )
            sleep(0.1)
            if printfig Print2Disk( f, path, "Vall", it) end
        end
    end
end

function Cerjean2D(X, Lbc, l, Δ)
    return ((1.0 .- exp.(-(X.x*ones(size(X.y))'.-0*l.x).^2/Lbc.^2))
         .*(1.0 .- exp.(-(X.x*ones(size(X.y))' .-  l.x).^2/Lbc.^2)) )
        #  .*(1.0 .- exp.(-(ones(size(X.x))*X.y' .-0*l.y).^2/Lbc.^2)) )
        #  .*(1.0 .- exp.(-(ones(size(X.x))*X.y' .-  l.y).^2/Lbc.^2)))
end

function Print2Disk( f, path, field, istep; res=4)
    path1 = path*"/_$field/"
    mkpath(path1)
    save(path1*"$field"*@sprintf("%05d", istep)*".png", f, px_per_unit = res) 
end

MainSource()
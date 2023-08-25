using SeismicQ, GLMakie, MathTeXEngine
Makie.update_theme!(fonts = (regular = texfont(), bold = texfont(:bold), italic = texfont(:italic)))

function MainRicker2D()
    # Spatial extent
    l  = (x = 25, y = 25)

    # Discretization
    Nc  = (x = 200, y = 200) 
    Δ   = (x = l.x/Nc.x, y = l.y/Nc.y, z=1.0)
    X   = (v = (x= LinRange(0,l.x,Nc.x+1)            , y= LinRange(0,l.y,Nc.y+1)),
            c = (x= LinRange(0-Δ.x/2,l.x+Δ.x/2,Nc.x+2) , y= LinRange(0-Δ.y/2,l.y+Δ.y/2,Nc.y+2)),
            i = (x= LinRange(0,l.x,Nc.x+1)             , y= LinRange(0-Δ.y/2,l.y+Δ.y/2,Nc.y+2)),
            j = (x= LinRange(0-Δ.x/2,l.x+Δ.x/2,Nc.x+2) , y= LinRange(0,l.y,Nc.y+1))) 
        
    # Array sizes
    szv   = (Nc.x+1, Nc.y+1)
    szc   = (Nc.x+2, Nc.y+2)

    # Source parameters
    𝑓₀   = 50   # Central frequency of the source [Hz]
    t₀   = 1.2/𝑓₀
    σ₀   = l.x/30
    x₀   = l.x/2
    y₀   = l.y/2
    t    = t₀    

    # Compute Ricker function with 2D spatial support
    f_ext = (v=zeros(szv)  , c=zeros(szc))
    x2d   = X.v.x * ones(size( X.v.y))'
    y2d   = ones(size( X.v.x)) * X.v.y'
    @. f_ext.v = Ricker.( x2d, x₀, y2d, y₀, t, t₀, 𝑓₀, σ₀)

    # Makie visualisation
    resol = 1000 
    f     = Figure(resolution = (l.x/l.y*resol, resol), fontsize=25)
    ax1   = Axis(f[1, 1], title = L" vz on v grid at $t$ = %$(t) [s]", xlabel = L"$x$ [m]", ylabel = L"$y$ [m]")
    hm    = GLMakie.heatmap!(ax1, X.v.x, X.v.y, f_ext.v, colormap = :turbo)
    colsize!(f.layout, 1, Aspect(1, l.x/l.y))
    GLMakie.Colorbar(f[1, 2], hm, label = "Vz [m/s]", width = 20, labelsize = 25, ticklabelsize = 14 )
    GLMakie.colgap!(f.layout, 20)
    # DataInspector(f)
    display(f)

end

MainRicker2D()
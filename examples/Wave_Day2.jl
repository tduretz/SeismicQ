using SeismicQ, Plots

function MainSource()
    
    # Spatial extent
    Lx  = 25.0

    # Mechanical parameters 
    ρ   = 1500.0
    K   = 1.e9
    G   = 1.e8
    


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
    # Courant criteria from wavespeed

    # Storage
    # ε̇xx
    # ∇v
    # P 
    # Vx
    # τxx
    time = zeros(Nt)
    acc  = zeros(Nt)
    vel  = zeros(Nt)
    v    = 0.

    # Time loop
    for it=1:Nt

        # Compute Ricker function
        t += Δt
        a  = Ricker(t, t₀, 𝑓₀)
        v += a*Δt

        # BC
        # ε̇xx
        # ∇v
        # τxx <------------ τxx = f(G, ε̇xx)
        # P   <------------   P = f(K,  ∇v)
        # Vx

        # For visualisation purpose
        time[it] = t
        acc[it]  = a
        vel[it]  = v
    end

    # Visualisation
    p1 = plot(time, acc, xlabel="t", ylabel="a")
    p2 = plot(time, vel, xlabel="t", ylabel="v")
    plot(p1, p2, layout=(2,1))

end

MainSource()
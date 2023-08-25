using SeismicQ, Plots,Classical

function MainSource()
    
    # Spatial extent
    L  = (x = 25.0, y = 12.5)

    # Mechanical parameters 
    ρ₀   = 1500.0
    K₀   = 1.e9
    G₀   = 1.e8
    c₀   = sqrt((K₀+4/3*G₀)/ρ₀) 
     
    # Discretization
    Nc  = (x = 100, y = 50) 
    Δ   = (x = L.x/Nc.x, y = L.y/Nc.y)
    X  = (v = (x= LinRange(0,L.x,Nc.x+1)            , y= LinRange(0,L.y,Nc.y+1)),
          c = (x= LinRange(0-Δ.x/2,L.x+Δ.x/2,Nc.x+2) , y= LinRange(0-Δ.y/2,L.y+Δ.y/2,Nc.y+2)),
          i = (x= LinRange(0,L.x,Nc.x+1)            , y= LinRange(0-Δ.y/2,L.y+Δ.y/2,Nc.y+2)),
          j = (x= LinRange(0-Δ.x/2,L.x+Δ.x/2,Nc.x+2) , y= LinRange(0,L.y,Nc.y+1))) 
    


    # Source parameters
    𝑓₀   = 50     # Central frequency of the source [Hz]
    t₀   = 1.2/𝑓₀
    isrc = Int((Nc.x/2)+1)
    jsrc = Int((Nc.y/2)+1)

    # Time domain
    Δt   = min(1e10, 0.3*Δ.x/c₀, 0.3*Δ.y/c₀ ) # Courant criteria from wavespeed
    Nt   = 2000
    Nout = 100
    t    = -t₀
   
    # Storage on centers # +2 for ghost nodes for BCs
    szv   = (Nc.x+1, Nc.y+1)
    szc   = (Nc.x+2, Nc.y+2)
    szi   = (Nc.x+1, Nc.y+2)
    szj   = (Nc.x+2, Nc.y+1)
    # Storage on i and j meshes
    K     = (i= ones(szi)*K₀, j= ones(szj)*K₀) 
    G     = (i= ones(szi)*G₀, j= ones(szj)*G₀) 
    ∇V    = (i = zeros(szi),  j = zeros(szj))
    P     = (i = zeros(szi),  j = zeros(szj))
    L     = (i = (xx=zeros(szi), xy=zeros(szi), yx=zeros(szi), yy=zeros(szi)),
             j = (xx=zeros(szj), xy=zeros(szj), yx=zeros(szj), yy=zeros(szj)))
             
    ε̇     = ( i=(xx=zeros(szi), yy=zeros(szi), zz=zeros(szi), xy=zeros(szi)),
              j=(xx=zeros(szj), yy=zeros(szj), zz=zeros(szj), xy=zeros(szj))) 
    
    τ     = ( i=(xx=zeros(szi), yy=zeros(szi), zz=zeros(szi), xy=zeros(szi)),
              j=(xx=zeros(szj), yy=zeros(szj), zz=zeros(szj), xy=zeros(szj)))           

    # Storage on v and c meshes
    V     = ( v=(x=zeros(szv), y=zeros(szv), z=zeros(szv)),
              c=(x=zeros(szc), y=zeros(szc), z=zeros(szc)))

    ρ     = (v=ones(szv)*ρ₀, c=ones(szc)*ρ₀)
    f_ext = (v=zeros(szv)  , c=zeros(szc))
    Vnorm = zeros(szc)
    # BC
    # Lbc        = 2.
    # # BC on v and c mesh
    # bc_filtW_V = (v= (1.0 .- exp.(-(X.v.x.-0L.x).^2/Lbc.^2)
    # bc_filtW_v = 1.0 .- exp.(-(X.v.x.-0L.x).^2/Lbc.^2)
    # bc_filtE_v = 1.0 .- exp.(-(xv.- L.x).^2/Lbc.^2)    
    # # BC on i and j mesh
    # bc_filtE_c = 1.0 .- exp.(-(xc.- Lx).^2/Lbc.^2)
    # bc_filtW_c = 1.0 .- exp.(-(xc.-0Lx).^2/Lbc.^2)

    # # Time loop
     @time for it=1:Nt

        # Compute Ricker function
        t          += Δt
        a           = Ricker(t, t₀, 𝑓₀)
        f_ext.v[isrc,jsrc] = ρ.v[isrc,jsrc]*a

        # Velocity gradient components
        #@show size(V.v.y[2:end,:])
        #@show size(L.j.x[2:end-1,2:end-1])

        @. L.i.xx[:,2:end-1] = (V.c.x[2:end,2:end-1] - V.c.x[1:end-1,2:end-1])/Δ.x
        @. L.j.xx[2:end-1,:] = (V.v.x[2:end,:] - V.v.x[1:end-1,:])/Δ.x

        @. L.i.yx[:,2:end-1] = (V.c.y[2:end,2:end-1] - V.c.y[1:end-1,2:end-1])/Δ.x
        @. L.j.yx[2:end-1,:] = (V.v.y[2:end,:] - V.v.y[1:end-1,:])/Δ.x

        @. L.i.yy[:,2:end-1] = (V.v.y[:,2:end] - V.v.y[:,1:end-1])/Δ.y
        @. L.j.yy[2:end-1,:] = (V.c.y[2:end-1,2:end] - V.c.y[2:end-1,1:end-1])/Δ.y

        @. L.i.xy[:,2:end-1] = (V.v.x[:,2:end] - V.v.x[:,1:end-1])/Δ.y
        @. L.j.xy[2:end-1,:] = (V.c.x[2:end-1,2:end] - V.c.x[2:end-1,1:end-1])/Δ.y
        
        
    #     # Divergence
        @. ∇V.i   = L.i.xx + L.i.yy
        @. ∇V.j   = L.j.xx + L.j.yy

    #     # Deviatoric strain rate 
        @. ε̇.i.xx = L.i.xx - 1//3*∇V.i
        @. ε̇.j.xx = L.j.xx - 1//3*∇V.j

        @. ε̇.i.yy = L.i.yy - 1//3*∇V.i
        @. ε̇.j.yy = L.j.yy - 1//3*∇V.j

        @. ε̇.i.zz = - 1//3*∇V.i
        @. ε̇.j.zz = - 1//3*∇V.j

        @. ε̇.i.xy = 1//2*(L.i.xy + L.i.yx)
        @. ε̇.j.xy = 1//2*(L.j.xy + L.j.yx)

      
    #     # Stress update
        @. τ.i.xx = f_shear(G.i)*Δt*(ε̇.i.xx) + f_relax(G.i)*τ.i.xx
        @. τ.j.xx = f_shear(G.j)*Δt*(ε̇.j.xx) + f_relax(G.j)*τ.j.xx

        @. τ.i.yy = f_shear(G.i)*Δt*(ε̇.i.yy) + f_relax(G.i)*τ.i.yy
        @. τ.j.yy = f_shear(G.j)*Δt*(ε̇.j.yy) + f_relax(G.j)*τ.j.yy

        @. τ.i.zz = f_shear(G.i)*Δt*(ε̇.i.zz) + f_relax(G.i)*τ.i.zz
        @. τ.j.zz = f_shear(G.j)*Δt*(ε̇.j.zz) + f_relax(G.j)*τ.j.zz

        @. τ.i.xy = f_shear(G.i)*Δt*(ε̇.i.xy) + f_relax(G.i)*τ.i.xy
        @. τ.j.xy = f_shear(G.j)*Δt*(ε̇.j.xy) + f_relax(G.j)*τ.j.xy

    #     # Pressure update 
        @. P.i    = P.i - Δt*f_bulk(K.i)*∇V.i
        @. P.j    = P.j - Δt*f_bulk(K.j)*∇V.j

    #     # Linear momentum balance
        @. V.v.x[2:end-1,2:end-1] = (V.v.x[2:end-1,2:end-1] 
                                    + Δt/ρ.v[2:end-1,2:end-1]
                                    *((τ.j.xx[3:end-1,2:end-1]-τ.j.xx[2:end-2,2:end-1])/Δ.x
                                    + (τ.i.xy[2:end-1,3:end-1]-τ.i.xy[2:end-1,2:end-2])/Δ.y 
                                    - (P.j[3:end-1,2:end-1]-P.j[2:end-2,2:end-1])/Δ.x 
                                    - f_ext.v[2:end-1,2:end-1]))
        
        @. V.c.x[2:end-1,2:end-1] = (V.c.x[2:end-1,2:end-1] 
                                    + Δt/ρ.c[2:end-1,2:end-1]
                                    *((τ.i.xx[2:end,2:end-1]-τ.i.xx[1:end-1,2:end-1])/Δ.x
                                    + (τ.j.xy[2:end-1,2:end]-τ.j.xy[2:end-1,1:end-1])/Δ.y 
                                    - (P.i[2:end,2:end-1]-P.i[1:end-1,2:end-1])/Δ.x 
                                    - f_ext.c[2:end-1,2:end-1]))                            

       @. V.v.y[2:end-1,2:end-1] = (V.v.y[2:end-1,2:end-1] 
                                    + Δt/ρ.v[2:end-1,2:end-1]
                                    *((τ.j.xy[3:end-1,2:end-1]-τ.j.xy[2:end-2,2:end-1])/Δ.x
                                    + (τ.i.yy[2:end-1,3:end-1]-τ.i.yy[2:end-1,2:end-2])/Δ.y 
                                    - (P.i[2:end-1,3:end-1]-P.i[2:end-1,2:end-2])/Δ.y 
                                    - f_ext.v[2:end-1,2:end-1]))
        
        @. V.c.y[2:end-1,2:end-1] = (V.c.y[2:end-1,2:end-1] 
                                    + Δt/ρ.c[2:end-1,2:end-1]
                                    *((τ.i.xy[2:end,2:end-1]-τ.i.xy[1:end-1,2:end-1])/Δ.x
                                    + (τ.j.yy[2:end-1,2:end]-τ.j.yy[2:end-1,1:end-1])/Δ.y 
                                    - (P.j[2:end-1,2:end]-P.j[2:end-1,1:end-1])/Δ.y 
                                    - f_ext.c[2:end-1,2:end-1]))   
       
    
    #     # Absorbing boundary Cerjean et al. (1985)
    #     @.  V.x  = V.x  * bc_filtW_v 
    #     @.  P    = P    * bc_filtW_c 
    #     @.  τ.xx = τ.xx * bc_filtW_c 
    #     @.  V.x  = V.x  * bc_filtE_v 
    #     @.  P    = P    * bc_filtE_c 
    #     @.  τ.xx = τ.xx * bc_filtE_c 

        # Visualisation
        if mod(it, Nout)==0
            @. Vnorm = sqrt(V.c.x^2+V.c.y^2)
            display(heatmap(X.c.x,X.c.y, Vnorm' ))
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
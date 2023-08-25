using SeismicQ, Plots

function MainSource()
    visu=true
    # Spatial extent
    l  = (x = 25.0, y = 12.5)

    # Mechanical parameters 
    ρ₀   = 1500.0
    K₀   = 1.e9
    G₀   = 1.e8
    c₀   = sqrt((K₀+4/3*G₀)/ρ₀) 
     
    # Discretization
    Nc  = (x = 100, y = 50) 
    Δ   = (x = l.x/Nc.x, y = l.y/Nc.y, z=1.0)
    X   = (v = (x= LinRange(0,l.x,Nc.x+1)            , y= LinRange(0,l.y,Nc.y+1)),
          c = (x= LinRange(0-Δ.x/2,l.x+Δ.x/2,Nc.x+2) , y= LinRange(0-Δ.y/2,l.y+Δ.y/2,Nc.y+2)),
          i = (x= LinRange(0,l.x,Nc.x+1)            , y= LinRange(0-Δ.y/2,l.y+Δ.y/2,Nc.y+2)),
          j = (x= LinRange(0-Δ.x/2,l.x+Δ.x/2,Nc.x+2) , y= LinRange(0,l.y,Nc.y+1))) 
    


    # Source parameters
    𝑓₀   = 100    # Central frequency of the source [Hz]
    t₀   = 1.2/𝑓₀
    src = (i=Int((Nc.x/2)+1),j=Int((Nc.y/2)+1))
    # src = (i=[Int(10/Δx) ],j=Int((Nc.y/2)+1))
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
    L     = (i = (xx=zeros(szi), xy=zeros(szi), yx=zeros(szi), yy=zeros(szi),zx=zeros(szi),zy=zeros(szi)),
             j = (xx=zeros(szj), xy=zeros(szj), yx=zeros(szj), yy=zeros(szj),zx=zeros(szj),zy=zeros(szj)))
             
    ε̇     = ( i=(xx=zeros(szi), yy=zeros(szi), zz=zeros(szi), xy=zeros(szi),xz=zeros(szi),yz=zeros(szi)),
              j=(xx=zeros(szj), yy=zeros(szj), zz=zeros(szj), xy=zeros(szj),xz=zeros(szj),yz=zeros(szj))) 
    
    τ     = ( i=(xx=zeros(szi), yy=zeros(szi), zz=zeros(szi), xy=zeros(szi),xz=zeros(szi),yz=zeros(szi)),
              j=(xx=zeros(szj), yy=zeros(szj), zz=zeros(szj), xy=zeros(szj),xz=zeros(szj),yz=zeros(szj)))           

    # Storage on v and c meshes
    V     = ( v=(x=zeros(szv), y=zeros(szv), z=zeros(szv)),
              c=(x=zeros(szc), y=zeros(szc), z=zeros(szc)))

    ρ     = (v=ones(szv)*ρ₀, c=ones(szc)*ρ₀)
    f_ext = (v=zeros(szv)  , c=zeros(szc))
    Vnorm = zeros(szc)
    # BC
     Lbc        = 1.
    # # BC on v and c mesh
     bc_filt_V   = (v=Cerjean2D(X.v,Lbc,l,Δ),c=Cerjean2D(X.c,Lbc,l,Δ))
     bc_filt_tau = (i=Cerjean2D(X.i,Lbc,l,Δ),j=Cerjean2D(X.j,Lbc,l,Δ))
     Vmax = 0.0
    #  bc_filtW_V   = (v=1.0 .- exp.(-(X.v.x*ones(size(X.v.y))' .-0*l.x).^2/Lbc.^2),
    #                  c=1.0 .- exp.(-(X.c.x*ones(size(X.c.y))' .-0*l.x).^2/Lbc.^2))
    #  bc_filtW_tau = (i=1.0 .- exp.(-(X.i.x*ones(size(X.i.y))' .-0l.x).^2/Lbc.^2),
    #                  j=1.0 .- exp.(-(X.j.x*ones(size(X.i.y))' .-0l.x).^2/Lbc.^2))
    # bc_filtE_v = 1.0 .- exp.(-(xv.- L.x).^2/Lbc.^2)    
    # # BC on i and j mesh
    # bc_filtE_c = 1.0 .- exp.(-(xc.- Lx).^2/Lbc.^2)
    # bc_filtW_c = 1.0 .- exp.(-(xc.-0Lx).^2/Lbc.^2)

    # # Time loop
     @views @time for it=1:Nt

        # Compute Ricker function
        t                  += Δt
        a                  = Ricker(t, t₀, 𝑓₀)
        # for isrc = 1:nsrc
        #    f_ext.v[src.i[isrc],src.j[isrc]] += ρ.v[src.i[isrc],src.j[isrc]]*a
        # end
        f_ext.v[src.i,src.j] = ρ.v[src.i,src.j]*a
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

        @. L.i.zy[:,2:end-1] = (V.v.z[:,2:end] - V.v.z[:,1:end-1])/Δ.y
        @. L.j.zy[2:end-1,:] = (V.c.z[2:end-1,2:end] - V.c.z[2:end-1,1:end-1])/Δ.y

        @. L.i.zx[:,2:end-1] = (V.c.z[2:end,2:end-1] - V.c.z[1:end-1,2:end-1])/Δ.x
        @. L.j.zx[2:end-1,:] = (V.v.z[2:end,:] - V.v.z[1:end-1,:])/Δ.x
        
        
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
        
               # in 2D Lxz and Lyz are zero 
        @. ε̇.i.xz = 1//2*(L.i.zx)
        @. ε̇.j.xz = 1//2*(L.j.zx)

        @. ε̇.i.yz = 1//2*(L.i.zy)
        @. ε̇.j.yz = 1//2*(L.j.zy)
      
    #     # Stress update
        @. τ.i.xx = f_shear(G.i)*Δt*(ε̇.i.xx) + f_relax(G.i)*τ.i.xx
        @. τ.j.xx = f_shear(G.j)*Δt*(ε̇.j.xx) + f_relax(G.j)*τ.j.xx

        @. τ.i.yy = f_shear(G.i)*Δt*(ε̇.i.yy) + f_relax(G.i)*τ.i.yy
        @. τ.j.yy = f_shear(G.j)*Δt*(ε̇.j.yy) + f_relax(G.j)*τ.j.yy

        @. τ.i.zz = f_shear(G.i)*Δt*(ε̇.i.zz) + f_relax(G.i)*τ.i.zz
        @. τ.j.zz = f_shear(G.j)*Δt*(ε̇.j.zz) + f_relax(G.j)*τ.j.zz

        @. τ.i.xy = f_shear(G.i)*Δt*(ε̇.i.xy) + f_relax(G.i)*τ.i.xy
        @. τ.j.xy = f_shear(G.j)*Δt*(ε̇.j.xy) + f_relax(G.j)*τ.j.xy

        @. τ.i.xz = f_shear(G.i)*Δt*(ε̇.i.xz) + f_relax(G.i)*τ.i.xz
        @. τ.j.xz = f_shear(G.j)*Δt*(ε̇.j.xz) + f_relax(G.j)*τ.j.xz

        @. τ.i.yz = f_shear(G.i)*Δt*(ε̇.i.yz) + f_relax(G.i)*τ.i.yz
        @. τ.j.yz = f_shear(G.j)*Δt*(ε̇.j.yz) + f_relax(G.j)*τ.j.yz

    #     # Pressure update 
        @. P.i    = P.i - Δt*f_bulk(K.i)*∇V.i
        @. P.j    = P.j - Δt*f_bulk(K.j)*∇V.j

    #     # Linear momentum balance
        @. V.v.x[2:end-1,2:end-1] = (V.v.x[2:end-1,2:end-1] 
                                    + Δt/ρ.v[2:end-1,2:end-1]
                                    *((τ.j.xx[3:end-1,2:end-1]-τ.j.xx[2:end-2,2:end-1])/Δ.x
                                    + (τ.i.xy[2:end-1,3:end-1]-τ.i.xy[2:end-1,2:end-2])/Δ.y 
                                    - (P.j[3:end-1,2:end-1]-P.j[2:end-2,2:end-1])/Δ.x 
                                    - 0.0*f_ext.v[2:end-1,2:end-1]))
        
        @. V.c.x[2:end-1,2:end-1] = (V.c.x[2:end-1,2:end-1] 
                                    + Δt/ρ.c[2:end-1,2:end-1]
                                    *((τ.i.xx[2:end,2:end-1]-τ.i.xx[1:end-1,2:end-1])/Δ.x
                                    + (τ.j.xy[2:end-1,2:end]-τ.j.xy[2:end-1,1:end-1])/Δ.y
                                    - (P.i[2:end,2:end-1]-P.i[1:end-1,2:end-1])/Δ.x 
                                    - 0.0*f_ext.c[2:end-1,2:end-1]))                            

       @. V.v.y[2:end-1,2:end-1] = (V.v.y[2:end-1,2:end-1] 
                                    + Δt/ρ.v[2:end-1,2:end-1]
                                    *((τ.j.xy[3:end-1,2:end-1]-τ.j.xy[2:end-2,2:end-1])/Δ.x
                                    + (τ.i.yy[2:end-1,3:end-1]-τ.i.yy[2:end-1,2:end-2])/Δ.y 
                                    - (P.i[2:end-1,3:end-1]-P.i[2:end-1,2:end-2])/Δ.y 
                                    - 0.0*f_ext.v[2:end-1,2:end-1]))
        
        @. V.c.y[2:end-1,2:end-1] = (V.c.y[2:end-1,2:end-1] 
                                    + Δt/ρ.c[2:end-1,2:end-1]
                                    *((τ.i.xy[2:end,2:end-1]-τ.i.xy[1:end-1,2:end-1])/Δ.x
                                    + (τ.j.yy[2:end-1,2:end]-τ.j.yy[2:end-1,1:end-1])/Δ.y 
                                    - (P.j[2:end-1,2:end]-P.j[2:end-1,1:end-1])/Δ.y 
                                    - 0.0*f_ext.c[2:end-1,2:end-1]))   

# the two terms in dPdz and dtauzzdz  cancel in linear elastic case ... but i am not sure with other rheologies so I have leavec them 
        @. V.v.z[2:end-1,2:end-1] = (V.v.z[2:end-1,2:end-1] 
                                    + Δt/ρ.v[2:end-1,2:end-1]
                                    *((τ.j.xz[3:end-1,2:end-1]-τ.j.xz[2:end-2,2:end-1])/Δ.x
                                    + (τ.i.yz[2:end-1,3:end-1]-τ.i.yz[2:end-1,2:end-2])/Δ.y 
                                    - f_ext.v[2:end-1,2:end-1]))
        
        @. V.c.z[2:end-1,2:end-1] = (V.c.z[2:end-1,2:end-1] 
                                    + Δt/ρ.c[2:end-1,2:end-1]
                                    *((τ.i.xz[2:end,2:end-1]-τ.i.xz[1:end-1,2:end-1])/Δ.x
                                    + (τ.j.yz[2:end-1,2:end]-τ.j.yz[2:end-1,1:end-1])/Δ.y 
                                    - f_ext.c[2:end-1,2:end-1]))   
    
    #     # Absorbing boundary Cerjean et al. (1985)
        
        @.  V.v.x  = V.v.x  * bc_filt_V.v 
        @.  V.v.y  = V.v.y  * bc_filt_V.v 
        @.  V.v.z  = V.v.z  * bc_filt_V.v
        @.  V.c.x  = V.c.x  * bc_filt_V.c 
        @.  V.c.y  = V.c.y  * bc_filt_V.c 
        @.  V.c.z  = V.c.z  * bc_filt_V.c  
        # @.  P    = P    * bc_filtW_c 
        # @.  τ.xx = τ.xx * bc_filtW_c 
    #     @.  V.x  = V.x  * bc_filtE_v 
    #     @.  P    = P    * bc_filtE_c 
    #     @.  τ.xx = τ.xx * bc_filtE_c 

        # Visualisation
        if mod(it, Nout)==0 && visu==true
           # @. Vnorm = sqrt(V.c.x^2+V.c.y^2)
           # Vmax = max(Vmax, maximum(V.v.z))
            display( heatmap(X.v.x,X.v.y, V.v.z' , clim=(-1.e-4,1.e-4)))
            sleep(0.1)
        end
    end
    #@show Vmax
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

function Cerjean2D(X,Lbc,l,Δ)
    return ((1.0 .- exp.(-(X.x*ones(size(X.y))'.-0*l.x).^2/Lbc.^2))
         .*(1.0 .- exp.(-(X.x*ones(size(X.y))' .-  l.x).^2/Lbc.^2))
         .*(1.0 .- exp.(-(ones(size(X.x))*X.y' .-0*l.y).^2/Lbc.^2))
         .*(1.0 .- exp.(-(ones(size(X.x))*X.y' .-  l.y).^2/Lbc.^2)))

        #  (1.0 .- exp.(-(X.v.x*ones(size(X.v.y))' .-0*l.x).^2/Lbc.^2))
        #  .*(1.0 .- exp.(-(X.v.x*ones(size(X.v.y))' .-l.x).^2/Lbc.^2))
        #  .*(1.0 .- exp.(-(X.v.x*ones(size(X.v.y))' .-0*l.y).^2/Lbc.^2))
        #  .*(1.0 .- exp.(-(X.v.x*ones(size(X.v.y))' .-l.y).^2/Lbc.^2))
end

MainSource()
@doc raw"""
    f_bulk(K) 

Volumetric part of the elastic constitutive update: 
provides the bulk modulus K for use in 

```math    
Δp = -K ∇⋅v Δt
```

"""
function f_bulk(K) 
    return K
end

"""
    f_shear(G)

Deviatoric part of the elastic constitutive update: 
provides the shear modulus G for use in 

```math    
Δτ = 2G Δε
```

"""
function f_shear(G)
    return G
end
 
# """
#     *(x, y, z...)

# Multiplication operator. `x * y * z *...` calls this function with multiple
# arguments, i.e. `*(x, y, z...)`.
# """
function f_relax(G)
    return 1.
end


"""
    f_shear(G)

Deviatoric part of the elastic constitutive update: 
provides the shear modulus G for use in 

```math    
Δτ = 2G Δε
```

"""
function f_shear(G,η,Δt) # Maxwell visco-elastic
    #return 2*G
    De=G*Δt/η
     return 2*G/(1+2*De)
    #return 2*G/(1+De)
end
 
function f_relax(G,η,Δt)
    #return 1-2*G*Δt/η
    De=G*Δt/η
    return 1/(1+2*De)
    #return (1-2*G*Δt)/(1-De)
end
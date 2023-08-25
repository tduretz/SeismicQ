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
    DeN=η/(G*Δt) # numerical Deborah number
    return 2.0*G*DeN/(1.0+DeN)
end

 
function f_relax(G,η,Δt)
    DeN=η/(G*Δt) # numerical Deborah number
    return DeN/(1.0+DeN)
end
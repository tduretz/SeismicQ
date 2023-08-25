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
 
"""
    f_relax(G)

Relaxation part of the deviatoric elastic constitutive update: 
by definition equal to 1 (no relaxation)

"""
function f_relax(G)
    return 1.
end


"""
    f_shear(G,η,Δt)

Deviatoric part of the Maxwell visco-elastic constitutive update: 
used to compute

```math    
Δτ = 2 G (DeN / (1+DeN)) Δε + Δτ_{RELAX}
```
with η the shear viscosity, and DeN a numerical Deborah number, equal to 
```math    
DeN = η / (G Δt)
```
"""
function f_shear(G,η,Δt) # Maxwell visco-elastic
    DeN=η/(G*Δt) # numerical Deborah number
    return 2.0*G*DeN/(1.0+DeN)
end

"""
f_relax(G,η,Δt)

Relaxation part of the deviatoric Maxwell visco-elastic constitutive update: 
for use in the calculation of 
```math    
Δτ = (DeN / (1+DeN)) τ_{OLD} + Δτ_{CONSTITUTIVE}
```
with DeN a numerical Deborah number, equal to 
```math    
DeN = η / (G Δt)
```
and η the shear viscosity.
"""
function f_relax(G,η,Δt)
    DeN=η/(G*Δt) # numerical Deborah number
    return DeN/(1.0+DeN)
end

"""
f_visc(ηb)

compute effective bulk viscosity for a Kelvin visco-elastic constitutive update: 
here obviously it is for a linear update 
"""
function f_visc(ηb)
    return ηb
end

"""
f_bulk(K,ηb,Δt)

compute effective bulk modulus for a Kelvin visco-elastic constitutive update 
    using the Fatboy number

"""
function f_bulk(K,ηb,Δt)
    Fb_bnum = ηb/K/Δt # numerical bulk fatboy number  
    return K*(1+Fb_bnum)
end
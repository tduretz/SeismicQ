#---------------------------------------------------#

@doc """
    θs(G,Δt)

Deviatoric part of the elastic constitutive update: 
provides the shear modulus G for use in 

```math    
Δτ = 2G Δε
```

"""
function θs(G,Δt)
    return 2.0*G*Δt
end

@doc """
    θs(G,η,Δt)

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
function θs(G,η,Δt) # Maxwell visco-elastic
    DeN = η/(G*Δt) # numerical Deborah number
    return 2.0*G*DeN/(1.0+DeN)*Δt
end

#---------------------------------------------------#

@doc raw"""
    χs(G,Δt)

Relaxation part of the deviatoric elastic constitutive update: 
by definition equal to 1 (no relaxation)

"""
function χs(G,Δt)
    return 1.
end

@doc raw"""
    χs(G,η,Δt)

Relaxation part of the deviatoric Maxwell visco-elastic constitutive update: 
for use in the calculation of 
```math    
Δτ = (DeN / (1+DeN)) τ_\mathrm{OLD} + Δτ_\mathrm{CONSTITUTIVE}
```
with DeN a numerical Deborah number, equal to 
```math    
DeN = η / (G Δt)
```
and η the shear viscosity.
"""
function χs(G,η,Δt)
    DeN = η/(G*Δt) # numerical Deborah number
    return DeN/(1.0+DeN)
end

#---------------------------------------------------#
#---------------------------------------------------#
#---------------------------------------------------#

@doc raw"""
    θb(K,Δt) 

Volumetric part of the elastic constitutive update: 
provides the bulk modulus K for use in 

```math    
Δp = -K ∇⋅v Δt
```

"""
function θb(K,Δt) 
    return -K*Δt
end

@doc raw"""
θb(K,ηb,Δt)

compute effective bulk viscosity for a Kelvin visco-elastic constitutive update 

"""
function θb(K,ηb,Δt)
    return -(Δt*K+ηb)
end

#---------------------------------------------------#

@doc raw"""
    χb(K,Δt)

compute effective bulk viscosity linear elastic update: 0.
"""
function χb(K,Δt) 
    return 0.
end

@doc raw"""
    χb(K,ηb,Δt)

compute effective bulk viscosity for a Kelvin visco-elastic constitutive update: 
here obviously it is for a linear update 
"""
function χb(K,ηb,Δt)
    return ηb
end

#---------------------------------------------------#
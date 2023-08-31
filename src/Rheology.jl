#---------------------------------------------------#

@doc """
    ηs(G,Δt)
    provides the effective shear viscosity ηs to model an elastic body. 
    In that case the of generic viscous constitutive update:
    ```math    
    τ =  θs(τ0-χs*ε̇0)+ 2 ηs ε̇   
    ``` 
    becomes 
    ```math    
    τ =  τ0 + 2 G Δt ε̇   
    ``` 
"""
function ηs(G,Δt)
    return 2.0*G*Δt
end

@doc """
    ηs(G,η,Δt)
    provides the effective shear viscosity ηs to model a Maxwell body. 
    In that case the of generic viscous constitutive update:
    ```math    
    τ =  θs(τ0-χs*ε̇0)+ 2 ηs ε̇   
    ``` 
    becomes 
    ```math    
    τ =  (DeN / (1+DeN)) τ0 + 2 G (DeN / (1+DeN))  Δt ε̇   
    ``` 
    with  DeN a numerical Deborah number, equal to 
    ```math    
    DeN = η / (G Δt)
    ```
    and 
    ```math
    ηs = 2 G (DeN / (1+DeN))  Δt
    ```
"""
function ηs(G,η,Δt) # Maxwell visco-elastic
    DeN = η/(G*Δt) # numerical Deborah number
    return 2.0*G*DeN/(1.0+DeN)*Δt
end

#---------------------------------------------------#

@doc raw"""
    χs(G,Δt) 
    In that case the of generic viscous constitutive update:
    ```math    
    τ =  θs(τ0-χs*ε̇0)+ 2 ηs ε̇   
    ``` 
    for elastic it becomes 
    ```math    
    τ =  τ0 + 2 G Δt ε̇   
    ``` 
    and 
    ```math    
    χs=0  
    ``` 
"""
function χs(G,Δt)
    return 0.
end

@doc raw"""
     χs(G,ηₘ,Δt) 
    In that case the of generic viscous constitutive update:
    ```math    
    τ =  θs(τ0-χs*ε̇0)+ 2 ηs ε̇   
    ``` 
    for maxwell elastic body it becomes 
    ```math    
    τ =  (DeN / (1+DeN)) τ0 + 2 G (DeN / (1+DeN))  Δt ε̇   
    ``` 
    with  DeN a numerical Deborah number, equal to 
    ```math    
    DeN = η / (G Δt)
    ```
    and 
    ```math    
    χs=0  
    ``` 
"""
function χs(G,ηₘ,Δt)
    return 0.
end
#---------------------------------------------------#
@doc raw"""
    θs(G,Δt)
    provides the Instaneous viscous part of τ to model an elastic body. 
    In that case the of generic viscous constitutive update:
    ```math   
    τ =  θs(τ0-χs*ε̇0)+ 2 ηs ε̇   
    ``` 
    becomes 
    ```math    
    τ =  τ0 + 2 G Δt ε̇   
    ``` 
    and 
    ```math    
    θs=1.  
    ``` 
"""
function θs(G,Δt)
    return G./G
end

@doc raw"""
    θs(G,η,Δt)
    provides the relaxation θs to model a Maxwell body. 
    In that case the of generic viscous constitutive update:
    ```math    
    τ =  θs(τ0-χs*ε̇0)+ 2 ηs ε̇   
    ``` 
    becomes 
    ```math    
    τ =  (DeN / (1+DeN)) τ0 + 2 G (DeN / (1+DeN))  Δt ε̇   
    ``` 
    with  DeN a numerical Deborah number, equal to 
    ```math    
    DeN = η / (G Δt)
    ```
    and 
    ```math
    θs = (DeN / (1+DeN))  
    ```
"""
function θs(G,η,Δt)
    DeN = η/(G*Δt) # numerical Deborah number
    return DeN/(1.0+DeN)
end

#---------------------------------------------------#
#------------------Bulk rheologies------------------#
#---------------------------------------------------#

@doc raw"""
    ηb(K,Δt) 
    Generic pressure update is  
    ```math    
    p = θb(p - (-χb∇⋅v0)) - ηb ∇⋅v 
    ```
    elastic pressure update is 
    ```math    
    p = p - K Δt ∇⋅v 
    ```
    with  
    ```math    
    ηb = K Δt
    ```
"""
function ηb(K,Δt) 
    return K*Δt
end

@doc raw"""
    ηb(K,ηₖ,Δt)
    Generic pressure update is  
    ```math    
    p = θb(p - (-χb∇⋅v0)) - ηb ∇⋅v 
    ```
    Kelvin visco-elastic pressure update is 
    ```math    
    p = (p - ηₖ ∇⋅v0) - (K Δt+ηₖ) ∇⋅v 
    ```
    with  
    ```math    
    ηb = K Δt+ηₖ
    ```
"""
function ηb(K,ηₖ,Δt)
    return Δt*K+ηₖ
end
#---------------------------------------------------#
@doc raw"""
    χb(K,Δt)
    Generic pressure update is  
    ```math    
    p = θb(p - (-χb∇⋅v0)) - ηb ∇⋅v 
    ```
    pressure update is 
    ```math    
    p = p - (K Δt) ∇⋅v 
    ```
    with  
    ```math    
    χb = 0.
    ```
    """
function χb(K,Δt) 
    return 0.
end

@doc raw"""
    χb(K,ηₖ,Δt)
    Generic pressure update is  
    ```math    
    p = θb(p - (-χb∇⋅v0)) - ηb ∇⋅v 
    ```
    Kelvin visco-elastic pressure update is 
    ```math    
    p = (p - ηₖ ∇⋅v0) - (K Δt+ηₖ) ∇⋅v 
    ```
    with  
    ```math    
    χb = ηₖ
    ```
"""
function χb(K,ηₖ,Δt)
    return ηₖ
end
#---------------------------------------------------#
#---------------------------------------------------#

@doc raw"""
    θb(K,Δt)
    Generic pressure update is  
    ```math    
    p = θb(p - (-χb∇⋅v0)) - ηb ∇⋅v 
    ```
    pressure update is 
    ```math    
    p = p - (K Δt) ∇⋅v 
    ```
    with  
    ```math    
    θb = 1.
    ```
"""
function θb(K,Δt) 
    return K./K
end

@doc raw"""
    θb(K,ηₖ,Δt)

    Generic pressure update is  

    ```math    
    p = θb(p - (-χb∇⋅v0)) - ηb ∇⋅v 
    ```
    Kelvin visco-elastic pressure update is 
    ```math    
    p = (p - ηₖ ∇⋅v0) - (K Δt+ηₖ) ∇⋅v 
    ```
    with  
    ```math    
    θb = 1.0
    ```
"""
function θb(K,ηₖ,Δt)
    return 1.
end
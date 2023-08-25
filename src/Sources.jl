@doc raw"""
    Ricker(t, t₀, 𝑓₀)

Compute the Ricker function for t, t₀ and 𝑓₀.

```math    
f = (1.0 - 2.0 (\pi 𝑓₀ (t - t₀))^2) \exp(-π^2 𝑓₀^2 (t - t₀)^2)
```

# Examples
```julia-repl
julia>  Ricker(0.0, 0.0, 0.0)
1.0
```
"""
Ricker(t, t₀, 𝑓₀) = (1.0 - 2.0 * (π * 𝑓₀ * (t - t₀))^2) * exp(-π^2 * 𝑓₀^2 * (t - t₀)^2) # Anonymous function

#----------------------------------------------------------------#

@doc raw"""
    Ricker(x, x₀, t, t₀, 𝑓₀)

Compute the Ricker function with a dirac.

```math    
f = \delta(x₀) \exp{(-((x-x₀)^2)/2.0/σ₀^2)} (1.0 - 2.0 (\pi 𝑓₀ (t - t₀))^2) \exp(-π^2.0𝑓₀^2.0(t - t₀)^2)
```

# Examples
```julia-repl
julia> Ricker(2., 2., 0.0, 0.0, 0.0)
0.1353352832366127
```
"""
Ricker(x, x₀, t, t₀, 𝑓₀) = (x==x₀) * Ricker(t, t₀, 𝑓₀)

#----------------------------------------------------------------#

@doc raw"""
    Ricker(x, x₀, t, t₀, 𝑓₀, σ₀)

Compute the Ricker function with 1D spatial support.

```math    
f = \exp{(-((x-x₀)^2)/2.0/σ₀^2)} (1.0 - 2.0 (\pi 𝑓₀ (t - t₀))^2) \exp(-π^2.0𝑓₀^2.0(t - t₀)^2)
```

# Examples
```julia-repl
julia> Ricker(2., 0., 0.0, 0.0, 0.0, 1)
0.1353352832366127
```
"""
Ricker(x, x₀, t, t₀, 𝑓₀, σ₀) = exp(-((x-x₀)^2)/2.0/σ₀^2) * Ricker(t, t₀, 𝑓₀)

#----------------------------------------------------------------#

@doc raw"""
    Ricker(x, x₀, y, y₀, t, t₀, 𝑓₀, σ₀)

Compute the Ricker function with 2D spatial support.

```math    
f = \exp{(-((x-x₀)^2)/2.0/σ₀^2)} (1.0 - 2.0 (\pi 𝑓₀ (t - t₀))^2) \exp(-π^2.0𝑓₀^2.0(t - t₀)^2)
```

# Examples
```julia-repl
julia> Ricker(2., 0., 2., 0., 0.0, 0.0, 0.0, 1)
0.01831563888873418
```
"""
Ricker(x, x₀, y, y₀, t, t₀, 𝑓₀, σ₀) = exp(-( (x-x₀)^2.0 + (y-y₀)^2.0)/2.0/σ₀^2) * Ricker(t, t₀, 𝑓₀)

#----------------------------------------------------------------#
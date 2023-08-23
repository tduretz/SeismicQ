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
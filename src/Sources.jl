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

@doc raw"""
    Gaussian(x, x₀, y, y₀, σ₀)

Compute the Gaussian function with 2D spatial support.

```math    
f = \exp(-((x - x₀)^2 + (y - y₀)^2) / (2.0 * σ₀^2))
```

"""
Gaussian(x, x₀, y, y₀, σ₀) =
    exp(-((x - x₀)^2 + (y - y₀)^2) / (2.0 * σ₀^2))


@doc raw"""
    ExplosiveX(x, x₀, y, y₀, t, t₀, 𝑓₀, σ₀)

Compute y component of explosive source mexican hat of given frequency with 2D spatial support.

```math    
f = ((x - x₀) / σ₀^2) * Gaussian(x, x₀, y, y₀, σ₀) * Ricker(t, t₀, 𝑓₀)
```

"""
ExplosiveX(x, x₀, y, y₀, t, t₀, 𝑓₀, σ₀) =
    ((x - x₀) / σ₀^2) * Gaussian(x, x₀, y, y₀, σ₀) * Ricker(t, t₀, 𝑓₀)

@doc raw"""
    ExplosiveY(x, x₀, y, y₀, t, t₀, 𝑓₀, σ₀)

Compute y component of explosive source mexican hat of given frequency with 2D spatial support.

```math    
f = ((x - x₀) / σ₀^2) * Gaussian(x, x₀, y, y₀, σ₀) * Ricker(t, t₀, 𝑓₀)
```

"""
ExplosiveY(x, x₀, y, y₀, t, t₀, 𝑓₀, σ₀) =
    ((y - y₀) / σ₀^2) * Gaussian(x, x₀, y, y₀, σ₀) * Ricker(t, t₀, 𝑓₀)


@doc raw"""
RotationalY(x, x₀, y, y₀, t, t₀, 𝑓₀, σ₀, θ)

Compute the y-component of an oriented pure rotational source.

"""
RotationalY(x, x₀, y, y₀, t, t₀, 𝑓₀, σ₀, θ) = begin
dx = x - x₀
dy = y - y₀
xp = cos(θ) * dx + sin(θ) * dy
G = Gaussian(x, x₀, y, y₀, σ₀)
Rt = Ricker(t, t₀, 𝑓₀)

dGdx = -(dx / σ₀^2) * G
dψdx = cos(θ) * G + xp * dGdx

-dψdx * Rt

end

@doc raw"""
RotationalX(x, x₀, y, y₀, t, t₀, 𝑓₀, σ₀, θ)

Compute the x-component of an oriented pure rotational source.

The source is derived from the stream function

"""
RotationalX(x, x₀, y, y₀, t, t₀, 𝑓₀, σ₀, θ) = begin
dx = x - x₀
dy = y - y₀
xp = cos(θ) * dx + sin(θ) * dy
G = Gaussian(x, x₀, y, y₀, σ₀)
Rt = Ricker(t, t₀, 𝑓₀)

dGdy = -(dy / σ₀^2) * G
dψdy = sin(θ) * G + xp * dGdy

dψdy * Rt

end

@doc raw"""
DoubleCoupleX(x, x₀, y, y₀, t, t₀, 𝑓₀, σ₀, θ)

Compute the x-component of a rotated double-couple source.
"""
DoubleCoupleX(x, x₀, y, y₀, t, t₀, 𝑓₀, σ₀, θ) = begin
G = Gaussian(x, x₀, y, y₀, σ₀)
dGdx = -((x - x₀) / σ₀^2) * G
dGdy = -((y - y₀) / σ₀^2) * G

Mxx = -sin(2θ)
Mxy =  cos(2θ)

(Mxx * dGdx + Mxy * dGdy) * Ricker(t, t₀, 𝑓₀)

end

@doc raw"""
DoubleCoupleY(x, x₀, y, y₀, t, t₀, 𝑓₀, σ₀, θ)

Compute the y-component of a rotated double-couple source.

"""
DoubleCoupleY(x, x₀, y, y₀, t, t₀, 𝑓₀, σ₀, θ) = begin
G = Gaussian(x, x₀, y, y₀, σ₀)
dGdx = -((x - x₀) / σ₀^2) * G
dGdy = -((y - y₀) / σ₀^2) * G

Myx =  cos(2θ)
Myy =  sin(2θ)

(Myx * dGdx + Myy * dGdy) * Ricker(t, t₀, 𝑓₀)

end
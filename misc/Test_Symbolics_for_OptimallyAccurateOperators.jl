# Sorry for the "bordel", Boris.
# The bunch of exercises (Wave_Day2_clean_super, Operators.jl, pascalTriangle, testOperators)
# are the contribution from Kensuke Wesly Inst. Konishi and me myself. 
# We are trying to seek high-oder operators for 1D,2D,3D differentiation adapted to FSG scheme here
# and apply a Geller-Takeuchi (1995,1998) scheme to find optimally accurate operators boundary
# evaluating the errors of mass matrix (acceleration term) and stiffness matrix (all the viso-elasto-whatever)
#
# Titi told me that I should learn how to code with Symbolics so let's go (it's 16h on Friday though)
#  
#  Nobu

#using Symbolics, SymbolicNumericIntegration
#using DifferentialEquations, ModelingToolkit

using Symbolics

#@parameters x₁ x₂ x₃
@variables ∂₁ ∂₂ ∂₃
@variables ΔX₁ ΔX₂ ΔX
@variables f
@variables m

# m-th order

@show 1//factorial(m)*∂₁^
for n=

@parameters t σ ρ β
@variables x(t) y(t) z(t)
D = Differential(t)

eqs = [D(D(x)) ~ σ * (y - x),
    D(y) ~ x * (ρ - z) - y,
    D(z) ~ x * y - β * z]

@named sys = ODESystem(eqs)
sys = structural_simplify(sys)

u0 = [D(x) => 2.0,
    x => 1.0,
    y => 0.0,
    z => 0.0]

p = [σ => 28.0,
    ρ => 10.0,
    β => 8 / 3]

tspan = (0.0, 100.0)
prob = ODEProblem(sys, u0, tspan, p, jac = true)
sol = solve(prob)
using Plots
plot(sol, idxs = (x, y))


# Symbols







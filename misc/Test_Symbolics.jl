using Symbolics, SymbolicNumericIntegration

# Symbols
@variables x, a, y

# Integrate
int = integrate(3x^3 + 2x - 5 + a, x)

# Show results
@show int[1]

# Generate function
f_expr = build_function(int[1], x, a )

# Evaluate function
f = eval(f_expr)

# Numerical application
f(1.0, 2.0)
import .Operators

n₁=8
n₂=8

Nₚ = (2*n₁+1)*(2*n₂+1)

# Order

numOrder = 8

numCoefs = (1+numOrder)*numOrder÷2


#𝑓=zeros(Float64,(2*n₁+1, 2*n₂+1))

Δx₁ = 1
Δx₂ = 1

PascalCoefs = zeros(Int64,(numOrder+1,numOrder))
NewCoefs = zeros(Rational,(numCoefs))

# Compute Pascal triangle


for i in 1:numOrder
    #@show i
    #@show Operators.pascal_triangle(i+1)
    
    #@show PascalCoefs[1:i+1,i]=Operators.pascal_triangle(i+1)
    #coefs[1:i,i] = Operators.pascal_triangle(i+1)
    #@show coefs[i,1:i]
end



@show Coefs=Operators.pascal_set(numOrder+1)


iCo=1
for i in 1:numOrder-1
    
    @show denominator=factorial(i)
    for j in 1:i+1
        @show iCo, Coefs[iCo], typeof(Coefs[iCo])
        NewCoefs[iCo]=Coefs[iCo] // denominator
        
        global iCo+=1
    end
end
@show NewCoefs



for ix₁=-n₁:n₁
    for ix₂=-n₂:-n₁

        Δx=ix₁*Δx₁
        Δy=ix₂*Δx₂

        index_X = n₁ + ix₁ + 1
        index_Y = n₂ + ix₂ + 1

        coefs[index_X,index_Y,1]=1
        coefs[index_X,index_Y,2]=Δx
        coefs[index_X,index_Y,3]=Δy
       





    end
end
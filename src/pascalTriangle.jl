function pascal_triangle(n)
    if n <= 0
        return [1]
    elseif n == 1
        return [1]
    else
        prev_row = pascal_triangle(n - 1)
        new_row = [1]
        for i in 1:(n - 2)
            push!(new_row, prev_row[i] + prev_row[i + 1])
        end
        push!(new_row, 1)
        return new_row
    end
end

function pascal_set(n)
    # pt = pascal_triangle(n)
    # @show pt
    v = pascal_triangle(1)
    #@show v
    for i in 2:n
        #@show i
        pt = pascal_triangle(i)
        # @show pt 
        v =(v, pt)
        v = vcat(v...)
        #@show v
    end
    @show v
end

function pascal_line(i, n, order = 5)
    pascals = pascal_set(order)
    #@show pascals
    pascal_size= size(pascals)[1]
    pascal_length =pascal_size * n 
    #@show pascal_length
    v = zeros(pascal_length)
    # pascal_set = pascal_set(i)
    for j in 1:pascal_size
        #@show j
        v[(i-1)*pascal_size + j] = pascals[j]
    end

    #@show v

end

function pascal_matrix(n, order = 5)
    pmat = transpose(pascal_line(1,n,order))
    # @show pmat 
    for i in 1:n
        pline = pascal_line(i,n,order)
        prow = transpose(pline)
        @show pmat

        pmat = [pmat; prow]
        # @show pline
        # @show prow

    end
    @show pmat
end



function print_pascal_triangle(n)
    for i in 1:n
        row = pascal_triangle(i)
        println(row)
    end
end


# print_pascal_triangle(10)

# pascal_set(5)

# pascal_line(3,3)
pascal_matrix(3)


# v1 = [1 2]
# v2 = [3 4]
# v = [v1;v2]
# @show v

# v = (v1, v2); 
# @show v
# v = hcat(v...)
# @show v

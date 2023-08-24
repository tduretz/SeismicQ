function pascal_triangle(n)
    if n <= 0
        return [1]
    elseif n == 1
        return [1, 1]
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

function print_pascal_triangle(n)
    for i in 1:n
        row = pascal_triangle(i)
        println(row)
    end
end

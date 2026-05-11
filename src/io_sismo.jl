using Printf

function make_filename(idf0, idDe)
    return @sprintf("f0-%03d_De-%03d.dat", idf0, idDe)
end

function write_sismo_file(filename, De, K0, G0, f0, rho0,
                           Lx, DeltaX, Nxs, Nt,
                           xsismo, time, velocity_matrix)

    open(filename, "w") do io

        println(io, "$De # De")
        println(io, "$K0 # K0")
        println(io, "$G0 # G0")
        println(io, "$f0 # f0")
        println(io, "$rho0 # rho0")
        println(io, "$Lx # Lx")
        println(io, "$DeltaX # DeltaX")
        println(io, "$Nxs # Nxs")
        println(io, "$Nt # Nt")

        println(io, "# BINARY_DATA_START")

        write(io, Float64.(xsismo))
        write(io, Float64.(time))
        write(io, Float64.(velocity_matrix))
    end
end

function read_sismo_file(filename)

    read_value(io, T) = begin
        line  = readline(io)
        value = split(line, "#")[1]
        parse(T, strip(value))
    end

    open(filename, "r") do io

        De      = read_value(io, Float64)
        K0      = read_value(io, Float64)
        G0      = read_value(io, Float64)
        f0      = read_value(io, Float64)
        rho0    = read_value(io, Float64)
        Lx      = read_value(io, Float64)
        DeltaX  = read_value(io, Float64)
        Nxs     = read_value(io, Int)
        Nt      = read_value(io, Int)

        marker = readline(io)
        @assert startswith(marker, "# BINARY_DATA_START")

        xsismo = Vector{Float64}(undef, Nxs)
        time   = Vector{Float64}(undef, Nt)
        velocity_matrix = Matrix{Float64}(undef, Nxs, Nt)

        read!(io, xsismo)
        read!(io, time)
        read!(io, velocity_matrix)

        return De, K0, G0, f0, rho0, Lx, DeltaX,
               Nxs, Nt,
               xsismo, time, velocity_matrix
    end
end
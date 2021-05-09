"""
    minishow(key, value)
Show an element of tuple as simple format.

# Examples
```julia
julia> A = (a = "1", b = 1, c = 1.0, d = -1.0)
julia> julia> Equate.minishow.(keys(A), values(A))
         a :          1 
         b :        1 
         c : 1.00000 
         d : -1.00000 
(nothing, nothing, nothing, nothing)
```
"""

function minishow(io::IO, k, v)
    if v isa Float64
        @printf io "%10s : %3.5f \n" k v
    elseif v isa Int64
        @printf io "%10s : %8i \n" k v
    elseif v isa String
        @printf io "%10s : %10s \n" k v
    else
        @printf io "%10s : can't be shown! \n" k
    end
end

function minishow(io::IO, A)
    minishow.(Ref(io), keys(A), values(A))
end

function Base.show(io::IO, x::FreqTab)
    println(io, "Frequency table stats.")
    minishow(io, x.stats)
    return
end

function Base.show(io::IO, x::NEATFreqTab)
    println(io, "Frequency table stats.")
    println(io, "* X (independent part)")
    minishow(io, x.statsX)
    println(io, "* V (common part)")
    minishow(io, x.statsV)
    return
end

function Base.show(io::IO, x::NEATEquateResult)
    println(io, "Equating design: NEAT\nEquated method: $(x.method).\nTo show the table, extract `table` element.")
    println(io, x.synthetic)
    return
end

function Base.show(io::IO, x::SGEquateResult)
    println(io, "Equating design: EG\nEquated method: $(x.method).\nTo show the table, extract `table` element.")
    return
end
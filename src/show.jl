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

function minishow(k, v)
    if v isa Float64
        @printf "%10s : %3.5f \n" k v
    elseif v isa Int64
        @printf "%10s : %8i \n" k v
    elseif v isa String
        @printf "%10s : %10s \n" k v
    else
        @printf "%10s : can't be shown! \n" k
    end
end

function minishow(A)
    minishow.(keys(A), values(A))
end

function Base.show(io::IO, x::FreqTab)
    println("Frequency table stats.")
    minishow(x.stats)
end

function Base.show(io::IO, x::NEATFreqTab)
    println("Frequency table stats.")
    println("* X (independent part)")
    minishow(x.statsX)
    println("* V (common part)")
    minishow(x.statsV)
end

function Base.show(io::IO, x::NEATEquateResult)
    println("Equating design: NEAT\nEquated method: $(x.method).\nTo show the table, extract `tbl` element.")
    println(x.synthetic)
end

function Base.show(io::IO, x::SGEquateResult)
    println("Equating design: EG\nEquated method: $(x.method).\nTo show the table, extract `tbl` element.")
end
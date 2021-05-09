using Bootstrap

## COPY 
function Base.copy(d::NamedTuple{(:X, :Y), Tuple{FreqTab, FreqTab}})
    return (X = copy(d.X), Y = copy(d.Y))
end

function Base.copy(t::FreqTab)
    FreqTab(
        copy(t.table),
        copy(t.raw), 
        copy(t.interval),
        merge(t.stats)
    )
end

function copy(x::SGEquateResult)
    SGEquateResult(
        x.method, 
        copy(x.table), 
        isnothing(x.estimates) ? nothing : merge(x.estimates), 
        copy(x.data)
    )
end

function Base.copy(d::NamedTuple{(:X, :Y), Tuple{NEATFreqTab, NEATFreqTab}})
    return (X = copy(d.X), Y = copy(d.Y))
end

function Base.copy(t::NEATFreqTab)
    NEATFreqTab(
        copy(t.tableX),
        copy(t.tableV),
        copy(t.rawX), 
        copy(t.rawV), 
        copy(t.intervalX),
        copy(t.intervalV),
        copy(t.marginal),
        merge(t.statsX),
        merge(t.statsV),
    )
end

function copy(x::NEATEquateResult)
    NEATEquateResult(
        x.method, 
        copy(x.table), 
        copy(x.synthetic), 
        isnothing(x.estimates) ? nothing : merge(x.estimates), 
        copy(x.data)
    )
end

##

function initialize!(t::FreqTab)
    t.table.freq .= zero(eltype(t.table.freq))
    t.table.cumfreq .= zero(eltype(t.table.cumfreq))
    t.table.prob .= zero(eltype(t.table.prob))
    t.table.cumprob .= zero(eltype(t.table.cumprob))
end

function initialize!(t::NEATFreqTab)
    # X
    t.tableX.freq .= zero(eltype(t.tableX.freq))
    t.tableX.cumfreq .= zero(eltype(t.tableX.cumfreq))
    t.tableX.prob .= zero(eltype(t.tableX.prob))
    t.tableX.cumprob .= zero(eltype(t.tableX.cumprob))
    # V
    t.tableV.freq .= zero(eltype(t.tableV.freq))
    t.tableV.cumfreq .= zero(eltype(t.tableV.cumfreq))
    t.tableV.prob .= zero(eltype(t.tableV.prob))
    t.tableV.cumprob .= zero(eltype(t.tableV.cumprob))
end



"""
    recalculate!(t::FreqTab)
Re-calculate `FreqTab` by using renewd raw vector. `interval` and the scale in table will not be changes.
"""
function recalculate!(t::FreqTab)
    t.stats = basicstats(t.raw)
    initialize!(t)
    freq = map(j -> count(i -> i == j, t.raw), t.table.scale)
    cumfreq = cumsum(freq)
    cumprob = cumsum(freq) ./ sum(freq)
    prob = freq ./ sum(freq)
    cumprob = cumprob
    for (i, s) in enumerate(freq)
        loc = [i]
        t.table[loc, :freq] .= freq[i]
        t.table[loc, :cumfreq] .= cumfreq[i]
        t.table[loc, :prob] .= prob[i]
        t.table[loc, :cumprob] .= cumprob[i]
    end
    return
end

function recalculate!(t::NEATFreqTab, form = ["X", "V"])
    t.statsX = basicstats(t.rawX)
    t.statsV = basicstats(t.rawV)
    initialize!(t)
    if "X" ∈ form
        # X
        freq = map(j -> count(i -> i == j, t.rawX), t.tableX.scale)
        cumfreq = cumsum(freq)
        cumprob = cumsum(freq) ./ sum(freq)
        prob = freq ./ sum(freq)
        cumprob = cumprob
        for (i, s) in enumerate(freq)
            loc = [i]
            t.tableX[loc, :freq] .= freq[i]
            t.tableX[loc, :cumfreq] .= cumfreq[i]
            t.tableX[loc, :prob] .= prob[i]
            t.tableX[loc, :cumprob] .= cumprob[i]
        end
    end
    if "V" ∈ form
        # V
        freq = map(j -> count(i -> i == j, t.rawV), t.tableV.scale)
        cumfreq = cumsum(freq)
        cumprob = cumsum(freq) ./ sum(freq)
        prob = freq ./ sum(freq)
        cumprob = cumprob
        for (i, s) in enumerate(freq)
            loc = [i]
            t.tableV[loc, :freq] .= freq[i]
            t.tableV[loc, :cumfreq] .= cumfreq[i]
            t.tableV[loc, :prob] .= prob[i]
            t.tableV[loc, :cumprob] .= cumprob[i]
        end
    end
    return 
end


# Attach methods to functions in  bootstrab.jl

function Bootstrap.draw!(x::NamedTuple{(:X, :Y), Tuple{T, T}}, o::NamedTuple{(:X, :Y), Tuple{T, T}}) where {T <: FreqTab}
    idx = sample(examineeID(x.X.raw), length(x.X.raw))
    idy = sample(examineeID(x.Y.raw), length(x.Y.raw))
    for (to, from) in enumerate(idx)
        o.X.raw[to] = x.X.raw[from]
    end
    for (to, from) in enumerate(idy)
        o.Y.raw[to] = x.Y.raw[from]
    end
    # Re-evaluate freqtab
    recalculate!(o.X)
    recalculate!(o.Y)
end

function Bootstrap.draw!(x::NamedTuple{(:X, :Y), Tuple{T, T}}, o::NamedTuple{(:X, :Y), Tuple{T, T}}) where {T <: NEATFreqTab}
    idx = sample(examineeID(x.X.rawX), length(x.X.rawX))
    idy = sample(examineeID(x.Y.rawX), length(x.Y.rawX))
    for (to, from) in enumerate(idx)
        o.X.rawX[to] = x.X.rawX[from]
        o.X.rawV[to] = x.X.rawV[from]
    end
    for (to, from) in enumerate(idy)
        o.Y.rawX[to] = x.Y.rawX[from]
        o.Y.rawV[to] = x.Y.rawV[from]
    end
    # Re-evaluate freqtab
    recalculate!(o.X)
    recalculate!(o.Y)
end

function Base.size(x::NamedTuple{(:X, :Y), Tuple{T, T}}) where {T <: FreqTab}
    "X $(x.X.stats.N)", "Y $(x.Y.stats.N)"
end

function Base.size(x::NamedTuple{(:X, :Y), Tuple{T, T}}) where {T <: NEATFreqTab}
    "X.Main $(x.X.statsX.N)", "X.Common $(x.X.statsV.N)", "Y.Main $(x.Y.statsX.N)", "Y.Main $(x.Y.statsV.N)"
end

# Standard Equating method



examineeID(X) = range(1, length = length(X))

"""
    bootSE()
Estimate bootstrap standard error of equating (SEE).

The definition of SEE is the standard deviation of equating results (i.e. In linear equating, the slope and intercept).
The source of SEE is random error in data. Furthermore, the source of random error is specified by the data collection design, the method of equatings and the sample size.
"""
function bootSE()

end
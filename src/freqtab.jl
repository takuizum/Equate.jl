# non equivalent common item design

"""
    FreqTab(tab, raw, interval stats)
The basic struct for all equating.

- `table` A data frame of test scales and thier frequaencies.
- `raw` A raw vector. The order of examinees should correspond to another `FreqTab.raw`
- `interval` The interval of the score scale.
- `stats` Basic stats of the frequency table.
"""
mutable struct FreqTab <: EG
    table
    raw
    interval
    stats
end

function basicstats(X)
    N = length(X)
    Nm = sum(X .=== missing)
    mins, maxs = minimum(X), maximum(X)
    (N = N, Missing = Nm, min = mins, maxs = maxs, μ = mean(X), σ = std(X), kurtosis = kurtosis(X), skewness = skewness(X))
end

# Frequaency table for equivalent group design
"""
    freqtab(X; interval = 1.0, scale = minimum(X):interval:maximum(X))
    
Create `FreqTab`, which is used for SG design methods in Equate package.

# Arguments

- `X` Vector of raw test score that dose not contain missing value. The order of examinees should correspond to another vector to be equated.
- `interval` The interval size of scale (must be Float64). Default is 1.0. If `scale` was specified by user, ignored.
- `scale` Vector or StepRange represents a scale of test score X.
"""
function freqtab(X; interval = 1.0, scale = minimum(X):interval:maximum(X))
    stats = basicstats(X)
    freq = map(j -> count(i -> i == j, X), scale)
    cumfreq = cumsum(freq)
    cumprob = cumsum(freq) ./ sum(freq)
    res = DataFrame(scale = scale, freq = freq, cumfreq = cumfreq, prob = freq ./ sum(freq), cumprob = cumprob)
    return FreqTab(res, X[X .!== missing], interval, stats)
end

# equivalent group design
mutable struct NEATFreqTab <: NEAT
    tableX
    tableV
    rawX # independent form
    rawV # common form
    intervalX
    intervalV
    marginal # conditional freqency
    statsX
    statsV
end
# Frequency table for nonequivalent gtoup design.
"""
    freqtab(X, V;intervalX = 1.0, intervalV = 1.0, scaleX = minimum(X):intervalX:maximum(X), scaleV = minimum(V):intervalV:maximum(V))

Create `NEATFreqTab` for NEAT design.

# Arguments

- `X` Vector of raw score, except missing values, of target test.
- `V` Vector of raw score, except missing values, of anchor test.

"""
function freqtab(X, V; intervalX = 1.0, intervalV = 1.0, scaleX = minimum(X):intervalX:maximum(X), scaleV = minimum(V):intervalV:maximum(V))
    if length(X) != length(V)
        println("X and V must be same length(test scores of which the same group).")
    end
    # missing data hangling
    if any(ismissing.(X)) || any(ismissing.(V))
        println("There are more than one missing value in X or V. Use the listwise deletion.")
        x_key = X .!== missing
        v_key = V .!== missing
        X = X[x_key .& x_key]
        V = V[v_key .& v_key]
    end
    statsX = basicstats(X)
    statsV = basicstats(V)
    freqx = map(j -> count(i -> i == j, X), scaleX)
    cumprobx = cumsum(freqx)./ sum(freqx)
    freqv = map(j -> count(i -> i == j, V), scaleV)
    cumprobv = cumsum(freqv)./ sum(freqv)
    tableX = DataFrame(scale = scaleX, freq = freqx, cumfreq = cumsum(freqx), cumprob = cumprobx, prob = freqx ./ sum(freqx))
    tableV = DataFrame(scale = scaleV, freq = freqv, cumfreq = cumsum(freqv), cumprob = cumprobv, prob = freqv ./ sum(freqv))
    marginaltable = zeros(Int64, length(scaleX), length(scaleV))
    for (xi, xv) in enumerate(scaleX), (vi, vv) in enumerate(scaleV)
        marginaltable[xi,vi] = count(i -> i == vv, V[X .== xv])
    end
    return NEATFreqTab(tableX, tableV, X, V, intervalX, intervalV, marginaltable, statsX, statsV)
end
# Frequency table from smoothed frequency table.
"""
    freqtab(X::EG, V::EG)

Create `NEATFreqTab` which has been smoothed (Log Linear or Kernel method).

# Arguments

- `X` Smoothed `FreqTab` of target test.
- `V` Smoothed `FreqTab` of anchor test.
"""
function freqtab(X::EG, V::EG)
    marginaltable = zeros(Int64, length(X.table.scale), length(V.table.scale))
    for (xi, xv) in enumerate(X.table.scale), (vi, vv) in enumerate(V.table.scale)
        marginaltable[xi,vi] = count(i -> i == vv, V.raw[X.raw .== xv])
    end
    return NEATFreqTab(X.table, V.table, X.raw, V.raw, X.interval, V.interval, marginaltable)
end

# non equivalent common item design
struct FreqTab <: EG
    tab::DataFrame
    raw::Vector
    interval::Float64
end
# Frequaency table for equivalent group design
function freqtab(X; interval = 1.0, scale = minimum(X):interval:maximum(X))
    freq = map(j -> count(i -> i == j, X), scale)
    cumfreq = cumsum(freq)
    cumprob = cumsum(freq) ./ sum(freq)
    res = DataFrame(scale = scale, freq = freq, cumfreq = cumfreq, prob = freq ./ sum(freq), cumprob = cumprob)
    return FreqTab(res, X, interval)
end

# equivalent group design
struct SGFreqTab <: NEAT
    tabX::DataFrame
    tabV::DataFrame
    rawX::Vector # independent form
    rawV::Vector # common form
    intervalX::Float64
    intervalV::Float64
    marginal::Matrix # conditional freqency
end
# Frequency table for nonequivalent gtoup design.
function freqtab(X, V;intervalX = 1.0, intervalV = 1.0, scaleX = minimum(X):intervalX:maximum(X), scaleV = minimum(V):intervalV:maximum(V))
    if length(X) != length(V)
        println("X and V must be same length(test scores of which the same group).")
    end
    freqx = map(j -> count(i -> i == j, X), scaleX)
    cumprobx = cumsum(freqx)./ sum(freqx)
    freqv = map(j -> count(i -> i == j, V), scaleV)
    cumprobv = cumsum(freqv)./ sum(freqv)
    tabX = DataFrame(scale = scaleX, freq = freqx, cumprob = cumprobx, prob = freqx ./ sum(freqx))
    tabV = DataFrame(scale = scaleV, freq = freqv, cumprob = cumprobv, prob = freqv ./ sum(freqv))
    marginaltable = zeros(Int64, length(scaleX), length(scaleV))
    for (xi, xv) in enumerate(scaleX), (vi, vv) in enumerate(scaleV)
        marginaltable[xi,vi] = count(i -> i == vv, V[X .== xv])
    end
    return SGFreqTab(tabX, tabV, X, V, intervalX, intervalV, marginaltable)
end
# Frequency table from smoothed frequency table.
function freqtab(X::EG, V::EG)
    marginaltable = zeros(Int64, length(X.tab.scale), length(V.tab.scale))
    for (xi, xv) in enumerate(X.tab.scale), (vi, vv) in enumerate(V.tab.scale)
        marginaltable[xi,vi] = count(i -> i == vv, V.raw[X.raw .== xv])
    end
    return SGFreqTab(X.tab, V.tab, X.raw, V.raw, X.interval, V.interval, marginaltable)
end

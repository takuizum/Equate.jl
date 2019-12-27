module Equate

using DataFrames, Statistics, GLM

export freqtab, round2, PRF, CDF, PFᵤ, PFₗ, e𝒀x, presmoothing

abstract type EquateDesign end
abstract type EGD <: EquateDesign end
abstract type NGD <: EquateDesign end
struct NGFreqTab <: EGD
    tab::DataFrame
    raw::Vector
    interval::Float64
end
struct EGFreqTab <: NGD
    tab::DataFrame
    rawX::Vector
    rawY::Vector
    intervalX::Float64
    intervalY::Float64
end
# Frequaency table for nonequivalent group design
function freqtab(X; interval = 1.0, scale = minimum(X):interval:maximum(X))
    freq = map(j -> count(i -> i == j, X), scale)
    cumfreq = cumsum(freq) ./ sum(freq)
    res = DataFrame(scale = scale, freq = freq, cumfreq = cumfreq)
    return NGFreqTab(res, X, interval)
end

# Frequency table for equivalent gtoup design
function freqtab(X, Y;intervalX = 1.0, intervalY = 1.0, scaleX = minimum(X):intervalX:maximum(X), scaleY = minimum(Y):intervalY:maximum(Y))
    freqx = map(j -> count(i -> i == j, X), scale)
    cumfreqx = cumsum(freqx)./ sum(freqx)
    freqy = map(j -> count(i -> i == j, Y), scale)
    cumfreqY = cumsum(freqy)./ sum(freqy)
    res = DataFrame(scale = scale, freqX = freqx, cumfreqX = cumfreqx, freqY = freqy, cumfreqY = cumfreqy)
    return EGFreqTab(res, X, Y, scaleX, scaleY, intervalX, intervalY)
end

# Equate method
abstract type EquateMethod end
struct Equipercentile <:EquateMethod end
struct Linear <: EquateMethod end
# Natural Round
round2(x; digits = 0) = sign(x) * floor( abs(x) * 10.0^digits + 0.5 ) / (10.0^digits)
# Percentile Rank Function
function CDF(x, F::NGFreqTab)
    if x < minimum(F.tab.scale) return 0 end
    if x > maximum(F.tab.scale) return 1 end
    F.tab.cumfreq[F.tab.scale .== x][1]
end
function PRF(x, F::NGFreqTab)
    if x < (minimum(F.tab.scale) - F.interval/2.0) return 0.0 end
    if x ≥ (maximum(F.tab.scale) + F.interval/2.0) return 100.0 end
    x⃰ = round2(x)
    Fx1 = CDF(x⃰-F.interval, F)#F.tab.cumfreq[F.tab.scale .== (x⃰-1.0)]
    Fx = CDF(x⃰, F)#F.tab.cumfreq[F.tab.scale .== x⃰]
    P = 100*(Fx1+(x-x⃰+F.interval/2.0)*(Fx-Fx1))[1]
    return P
end
# Percentile Function
function p_search_descend(P, F::NGFreqTab, offset)
    x = nothing;iter = length(F.tab.scale)
    while x == nothing
        iter -= 1
        x =  100CDF(F.tab.scale[iter], F) > P ? nothing : F.tab.scale[iter+offset]
    end
    return x
end
function p_search_ascend(P, F::NGFreqTab, offset)
    x = nothing;iter = 0
    while x == nothing
        iter += 1
        x = 100CDF(F.tab.scale[iter], F) < P ? nothing : iter == 1 ? 0.0 : F.tab.scale[iter+offset]
    end
    return x
end

function PFᵤ(P, F::NGFreqTab)
    if P ≥ 100.0 return (maximum(F.tab.scale) + .5) end
    xᵤ = P > 50.0 ? p_search_descend(P, F, 1) : p_search_ascend(P, F, 0)
    x = (P/100 - CDF(xᵤ-F.interval, F)) / (CDF(xᵤ, F) - CDF(xᵤ-F.interval, F))
    return isinf(x) || isnan(x) ? xᵤ -F.interval/2.0 : x + xᵤ -F.interval/2.0
end
function PFₗ(P, F::NGFreqTab)
    if P ≤ 0.0 return -.5 end
    xₗ = P > 50.0 ? p_search_descend(P, F, 0) : p_search_ascend(P, F, -1)
    x = (P/100 - CDF(xₗ, F)) / (CDF(xₗ+F.interval, F) - CDF(xₗ, F))
    return isinf(x) ? xₗ + F.interval/2.0 : x + xₗ + F.interval/2.0
end

function e𝒀x(X::NGFreqTab, Y::NGFreqTab; case = :upper)
    scaleY = Y.tab.scale
    eYxᵤ = zeros(Float64, length(scaleY)); eYxₗ = zeros(Float64, length(scaleY))
    for (i,v) in enumerate(scaleY)
        P = PRF(v, X)
        eYxᵤ[i] = PFᵤ(P, Y)
        eYxₗ[i] = PFₗ(P, Y)
    end
    if case == :upper
        eYx = eYxᵤ
    elseif case == :lower
        eYx = eYxₗ
    elseif case == :both
        eYx = string.(eYxᵤ, "_", eYxₗ)
    elseif case == :middle
        eYx = (eYxᵤ .+ eYxₗ) ./ 2.0
    end
    return DataFrame(scaleY = scaleY, eYx = eYx)
end

# Equating Function
function equate(X::EGFreqTab, method::Equipercentile)
    printle("hoge")
end

# LogLinear Transformation
function presmoothing(F::NGFreqTab; df = 4)
    # Creata GLM formula
    fml = "@formula(freq ~ "
    for d in 1:df
        fml *= "scale^$d"
        if d != df
            fml *= " + "
        else
            fml *= ")"
        end
    end
    println(fml)
    fit1 = glm(eval(Meta.parse(fml)), F.tab, Poisson(), LogLink())
    return fit
end

#-----------------
end # module

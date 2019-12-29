module Equate

using DataFrames, Statistics, GLM

export freqtab, round2, PRF, CDF, PFu, PFl, equipercentile, presmoothing, linear, Tuker, chainedlinear
export LogLinearFormula

abstract type EquateDesign end
abstract type EGD <: EquateDesign end
abstract type NGD <: EquateDesign end
# non equivalent common item design
struct FreqTab <: EGD
    tab::DataFrame
    raw::Vector
    interval::Float64
end
# Frequaency table for equivalent group design
function freqtab(X; interval = 1.0, scale = minimum(X):interval:maximum(X))
    freq = map(j -> count(i -> i == j, X), scale)
    cumfreq = cumsum(freq) ./ sum(freq)
    res = DataFrame(scale = scale, freq = freq, cumfreq = cumfreq)
    return FreqTab(res, X, interval)
end

# Equate method
abstract type EquateMethod end
struct Equipercentile <:EquateMethod end
struct Linear <: EquateMethod end
# Natural Round
round2(x; digits = 0) = sign(x) * floor( abs(x) * 10.0^digits + 0.5 ) / (10.0^digits)
# Percentile Rank Function
function CDF(x, F::FreqTab)
    if x < minimum(F.tab.scale) return 0 end
    if x > maximum(F.tab.scale) return 1 end
    F.tab.cumfreq[F.tab.scale .== x][1]
end
function PRF(x, F::FreqTab)
    if x < (minimum(F.tab.scale) - F.interval/2.0) return 0.0 end
    if x ≥ (maximum(F.tab.scale) + F.interval/2.0) return 100.0 end
    x′ = round2(x)
    Fx1 = CDF(x′-F.interval, F)#F.tab.cumfreq[F.tab.scale .== (x⃰-1.0)]
    Fx = CDF(x′, F)#F.tab.cumfreq[F.tab.scale .== x⃰]
    P = 100*(Fx1+(x-x′+F.interval/2.0)*(Fx-Fx1))[1]
    return P
end
# Percentile Function
function p_search_descend(P, F::FreqTab, offset)
    x = nothing;iter = length(F.tab.scale)
    while x == nothing
        iter -= 1
        x =  100CDF(F.tab.scale[iter], F) > P ? nothing : F.tab.scale[iter+offset]
    end
    return x
end
function p_search_ascend(P, F::FreqTab, offset)
    x = nothing;iter = 0
    while x == nothing
        iter += 1
        x = 100CDF(F.tab.scale[iter], F) < P ? nothing : iter == 1 ? 0.0 : F.tab.scale[iter+offset]
    end
    return x
end
function PFu(P, F::FreqTab)
    if P ≥ 100.0 return (maximum(F.tab.scale) + .5) end
    xu = P > 50.0 ? p_search_descend(P, F, 1) : p_search_ascend(P, F, 0)
    x = (P/100 - CDF(xu-F.interval, F)) / (CDF(xu, F) - CDF(xu-F.interval, F))
    return isinf(x) || isnan(x) ? xu -F.interval/2.0 : x + xu -F.interval/2.0
end
function PFl(P, F::FreqTab)
    if P ≤ 0.0 return -.5 end
    xl = P > 50.0 ? p_search_descend(P, F, 0) : p_search_ascend(P, F, -1)
    x = (P/100 - CDF(xl, F)) / (CDF(xl+F.interval, F) - CDF(xl, F))
    return isinf(x) ? xl + F.interval/2.0 : x + xl + F.interval/2.0
end
# equipercentile equating
function equipercentile(X::FreqTab, Y::FreqTab; case = :upper)
    scaleY = Y.tab.scale
    eYxu = zeros(Float64, length(scaleY)); eYxl = zeros(Float64, length(scaleY))
    for (i,v) in enumerate(scaleY)
        P = PRF(v, X)
        eYxu[i] = PFu(P, Y)
        eYxl[i] = PFl(P, Y)
    end
    if case == :upper
        eYx = eYxu
    elseif case == :lower
        eYx = eYxl
    elseif case == :both
        eYx = string.(eYxu, "_", eYxl)
    elseif case == :middle
        eYx = (eYxu .+ eYxl) ./ 2.0
    end
    return DataFrame(scaleY = scaleY, eYx = eYx)
end
# linear equating
function linear(X::FreqTab, Y::FreqTab)
    μX = mean(X.raw); σX = std(X.raw)
    μY = mean(Y.raw); σY = std(Y.raw)
    slope = σY/σX; intercept = μY - slope*μX
    eYx = @. X.tab.scale * slope + intercept
    return DataFrame(scaleX = X.tab.scale, eYx = eYx)
end
# LogLinear Transformation
function LogLinearFormula(df)
    fml = "@formula(freq ~ "
    for d in 1:df
        fml *= "scale^$d"
        if d != df
            fml *= " + "
        else
            fml *= ")"
        end
    end
    return eval(Meta.parse(fml))
end

function presmoothing(F::FreqTab; fml = LogLinearFormula(4))
    fit1 = glm(fml, F.tab, Poisson(), LogLink())
    pred = predict(fit1, DataFrame(scale = F.tab.scale))
    tab = DataFrame(scale = F.tab.scale, freq = pred, cumfreq = cumsum(pred))
    return FreqTab(tab, F.raw, F.interval), fit1
end
# equivalent group design
struct NEGFreqTab <: NGD
    tabX::DataFrame
    tabV::DataFrame
    rawX::Vector # independent form
    rawV::Vector # common form
    intervalX::Float64
    intervalV::Float64
    marginal::Matrix # conditional freqency
end
# Frequency table for nonequivalent gtoup design
function freqtab(X, V;intervalX = 1.0, intervalV = 1.0, scaleX = minimum(X):intervalX:maximum(X), scaleV = minimum(V):intervalV:maximum(V))
    if length(X) != length(V)
        println("X and V must be same length(test scores of which the same group).")
    end
    freqx = map(j -> count(i -> i == j, X), scaleX)
    cumfreqx = cumsum(freqx)./ sum(freqx)
    freqv = map(j -> count(i -> i == j, V), scaleV)
    cumfreqv = cumsum(freqv)./ sum(freqv)
    tabX = DataFrame(scale = scaleX, freq = freqx, cumfreq = cumfreqx)
    tabV = DataFrame(scale = scaleV, freq = freqv, cumfreq = cumfreqv)
    marginaltable = zeros(Int64, length(scaleX), length(scaleV))
    for (xi, xv) in enumerate(scaleX), (vi, vv) in enumerate(scaleV)
        marginaltable[xi,vi] = count(i -> i == vv, V[X .== xv])
    end
    return NEGFreqTab(tabX, tabV, X, V, intervalX, intervalV, marginaltable)
end
# Nonequivalent Groups : Linear methods
function Tuker(X::NEGFreqTab, Y::NEGFreqTab)
    # synsetic weight
    w₁ = length(X.raw) / (length(X.rawX) + length(X.rawV))
    w₂ = 1.0 - w1
    # test score
    x = X.rawX; xv = X.rawV
    y = Y.rawX; yv = Y.rawV
    # regression slope
    γ₁ = cov(x,xy) / var(xv)
    γ₂ = cov(y,yv) / var(yv)
    # synsetic mean and var
    μsX = mean(x) - w₂*gamma₁*(mean(xv)-mean(yv))
    μsY = mean(y) - w₁*gamma₂*(mean(xv)-mean(yv))
    σ²sX = var(x) - w₂*gamma₁^2*(var(xv)-var(yv)) + w₁*w₂*γ₁^2*(mean(xv)-mean(yv))^2
    σ²sY = var(y) - w₁*gamma₂^2*(var(xv)-var(yv)) + w₁*w₂*γ₂^2*(mean(xv)-mean(yv))^2
    # transformation
    slope = sqrt(σ²sY)/sqrt(σ²sX); intercept = μsY - slope*μsX
    eYx = @. X.tabX.scale * slope + intercept
    return DataFrame(scaleX = X.tabX.scale, eYx = eYx)
end
# Nonequivalent Goups : Chained linear Observed Score Equating
function chainedlinear(X::NEGFreqTab, Y::NEGFreqTab)
    # ******************************************** #
    # 1. put X on the scale of V -call lV(x);
    # 2. put V on the scale of Y - call lY(v);
    # 3. obtain Y-equivalent as lY(x).
    # ******************************************** #
    # test score
    x = X.rawX; xv = X.rawV
    y = Y.rawX; yv = Y.rawV
    # regression slope
    γ₁ = cov(x,xy) / var(xv)
    γ₂ = cov(y,yv) / var(yv)
    # moments
    μx = mean(x); σx = std(x)
    μxv = mean(xv); σxv = std(xv)
    μy = mean(y); σy = std(y)
    μyv = mean(yv); σyv = std(yv)
    slope = (σy/σyv)/(σx/σxv)
    intercept = μy + σy/σyv *(μxv - μyv) - slope * μx
    eYx = @. X.tabX.scale * slope + intercept
    return DataFrame(scaleX = X.tabX.scale, eYx = eYx)
end


#-----------------
end # module

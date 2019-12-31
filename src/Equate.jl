module Equate

using DataFrames, Statistics, GLM, Distributions, Optim

export freqtab, round2, PRF, CDF, PFu, PFl, Equipercentile, presmoothing, Linear, Tucker, ChainedLinear, ObservableStats, FrequencyEstimation, KernelSmoothing, BandwidthPenalty, EstBandwidth, BraunHolland, ChainedEquipercentile
export LogLinearFormula

abstract type EquateDesign end
abstract type EG <: EquateDesign end
abstract type NEAT <: EquateDesign end
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

# Equate method
abstract type SGEquateMethod end
# Natural Round
round2(x; digits = 0) = sign(x) * floor( abs(x) * 10.0^digits + 0.5 ) / (10.0^digits)
# Percentile Rank Function
function CDF(x, F::EG)
    if x < minimum(F.tab.scale) return 0 end
    if x > maximum(F.tab.scale) return 1 end
    F.tab.cumprob[F.tab.scale .== x][1]
end
function PRF(x, F::EG)
    if x < (minimum(F.tab.scale) - F.interval/2.0) return 0.0 end
    if x ≥ (maximum(F.tab.scale) + F.interval/2.0) return 100.0 end
    x′ = round2(x)
    Fx1 = CDF(x′-F.interval, F)#F.tab.cumfreq[F.tab.scale .== (x⃰-1.0)]
    Fx = CDF(x′, F)#F.tab.cumfreq[F.tab.scale .== x⃰]
    P = 100*(Fx1+(x-x′+F.interval/2.0)*(Fx-Fx1))[1]
    return P
end
# Percentile Function
function p_search_descend(P, F::EG, offset)
    x = nothing;iter = length(F.tab.scale)
    while x == nothing
        iter -= 1
        x =  100CDF(F.tab.scale[iter], F) > P ? nothing : F.tab.scale[iter+offset]
    end
    return x
end
function p_search_ascend(P, F::EG, offset)
    x = nothing;iter = 0
    while x == nothing
        iter += 1
        x = 100CDF(F.tab.scale[iter], F) < P ? nothing : iter == 1 ? 0.0 : F.tab.scale[iter+offset]
    end
    return x
end
function PFu(P, F::EG)
    if P ≥ 100.0 return (maximum(F.tab.scale) + .5) end
    xu = P > 50.0 ? p_search_descend(P, F, 1) : p_search_ascend(P, F, 0)
    x = (P/100 - CDF(xu-F.interval, F)) / (CDF(xu, F) - CDF(xu-F.interval, F))
    return isinf(x) || isnan(x) ? xu -F.interval/2.0 : x + xu -F.interval/2.0
end
function PFl(P, F::EG)
    if P ≤ 0.0 return -.5 end
    xl = P > 50.0 ? p_search_descend(P, F, 0) : p_search_ascend(P, F, -1)
    x = (P/100 - CDF(xl, F)) / (CDF(xl+F.interval, F) - CDF(xl, F))
    return isinf(x) ? xl + F.interval/2.0 : x + xl + F.interval/2.0
end
# equipercentile equating
struct ResultEquipercentile <: SGEquateMethod
    table::DataFrame
end
function Equipercentile(X::EG, Y::EG; case = :middle)
    scaleY = Y.tab.scale
    eYxu = zeros(Float64, length(scaleY)); eYxl = zeros(Float64, length(scaleY))
    for (i,v) in enumerate(scaleY)
        P = PRF(v, Y)
        eYxu[i] = PFu(P, X)
        eYxl[i] = PFl(P, X)
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
    tbl = DataFrame(scaleY = scaleY, eYx = eYx)
    return ResultEquipercentile(tbl)
end
# linear equating
struct ResultLinear <: SGEquateMethod
    table::DataFrame
end
function Linear(X::EG, Y::EG)
    μX = mean(X.raw); σX = std(X.raw)
    μY = mean(Y.raw); σY = std(Y.raw)
    slope = σY/σX; intercept = μY - slope*μX
    eYx = @. X.tab.scale * slope + intercept
    tbl = DataFrame(scaleX = X.tab.scale, eYx = eYx)
    ResultLinear(tbl)
end
# LogLinear Transformation
function LogLinearFormula(df)
    fml = "@formula(freq ~ 1 +"
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
struct SmoothedFreqTab <: EG
    tab::DataFrame
    raw::Vector
    interval::Float64
end
function presmoothing(F::FreqTab; fml = LogLinearFormula(4))
    fit1 = glm(fml, F.tab, Poisson(), LogLink())
    pred = predict(fit1, DataFrame(scale = F.tab.scale))
    tab = DataFrame(scale = F.tab.scale, prob = pred, cumprob = cumsum(pred))
    return SmoothedFreqTab(tab, F.raw, F.interval), fit1
end
# Kernel method
function RjX(x, xⱼ, a, μ, hX)
    return (x - a * xⱼ - (1-a)*μ) / (a*hX)
end
struct KernelFreqTab <: EG
    tab::DataFrame
    raw::Vector
    interval::Float64
    Bandwidth::Float64
end
function KernelSmoothing(X::EG; kernel = :Gaussian, hX = 0.66, scale = X.tab.scale)
    # hX = bandwidht of cumulative distribution function
    μ = mean(X.raw); σ² = var(X.raw)
    a² = σ² / (σ² + hX^2)
    a =sqrt(a²)
    𝒇hX = zeros(Float64, length(scale))
    FhX = zeros(Float64, length(scale))
    for (i, x) in enumerate(scale), (j, xⱼ) in enumerate(X.tab.scale)
        𝒇hX[i] += X.tab.prob[j]*pdf.(Normal(0, 1), RjX(x, xⱼ, a, μ, hX)) / (a*hX)
        FhX[i] += X.tab.prob[j]*cdf.(Normal(0, 1), RjX(x, xⱼ, a, μ, hX))
    end
    tbl = DataFrame(scale = scale, prob = 𝒇hX, cumprob = cumsum(𝒇hX))
    return KernelFreqTab(tbl, X.raw, X.interval, hX)
end
function BandwidthPenalty(hX, X::EG; kernel = :Gaussian)
    hX = exp(hX[1])
    r = X.tab.prob
    μ = mean(X.raw); σ² = var(X.raw)
    a² = σ² / (σ² + hX^2)
    a =sqrt(a²)
    𝒇hX = zeros(Float64, length(X.tab.scale))
    𝒇′hX = zeros(Float64, length(X.tab.scale))
    for (i, x) in enumerate(X.tab.scale), (j, xⱼ) in enumerate(X.tab.scale)
        R = RjX(x, xⱼ, a, μ, hX)
        𝒇hX[i] += X.tab.prob[j]*pdf.(Normal(0, 1), R) / (a*hX)
        𝒇′hX[i] -= X.tab.prob[j]*pdf.(Normal(0, 1), R) / (a*hX)^2 * R
    end
    pen1 = sum((r .- 𝒇hX) .^2)
    return pen1
end
function EstBandwidth(X::EG; kernel = :Gaussian)
    opt = optimize(hX -> BandwidthPenalty(hX, X), [0.5], method = BFGS())
    println("Minimizer sould be transformed `exp()` before interpretation. Minimizer = $(exp(opt.minimizer[1]))")
    return opt
end


# equivalent group design
abstract type NEATEquateMethod end
struct SGFreqTab <: NEAT
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
# Nonequivalent Groups : Linear methods
function ObservableStats(F::NEAT)
    x = F.rawX; v = F.rawV
    μx = mean(x); σx = std(x)
    μv = mean(v); σv = std(v)
    covxv = cov(x, v); corxv = cor(x, v)
    return μx, σx, μv, σv, covxv, corxv
end
struct ResultTucker <: NEATEquateMethod
    table::DataFrame
    synthetic::DataFrame
    estimates::NamedTuple
end
function Tucker(X::NEAT, Y::NEAT; w₁ = length(X.rawX) / (length(X.rawX) + length(Y.rawX)), w₂ = 1.0 - w₁)
    W = w₁ + w₂
    w₁ = w₁ / W; w₂ = w₂ / W
    # test score
    x = X.rawX; xv = X.rawV
    y = Y.rawX; yv = Y.rawV
    # ObservableStats
    μx, σx, μxv, σxv, covxv, corxv = ObservableStats(X)
    μy, σy, μyv, σyv, covyv, coryv = ObservableStats(Y)
    # regression slope
    γ₁ = covxv / σxv^2
    γ₂ = covyv / σyv^2
    # synthetic mean and var
    μsX = μx - w₂*γ₁*(μxv-μyv)
    μsY = μy + w₁*γ₂*(μxv-μyv)
    σ²sX = σx^2 - w₂*γ₁^2*(σxv^2-σyv^2) + w₁*w₂*γ₁^2*(μxv-μyv)^2
    σ²sY = σy^2 + w₁*γ₂^2*(σxv^2-σyv^2) + w₁*w₂*γ₂^2*(μxv-μyv)^2
    # transformation
    slope = sqrt(σ²sY)/sqrt(σ²sX); intercept = μsY - slope*μsX
    lYx = @. X.tabX.scale * slope + intercept
    tbl = DataFrame(scaleX = X.tabX.scale, lYx = lYx)
    return ResultTucker(tbl,
                       DataFrame(Group = [1, 2], μ = [μsX, μsY], σ = [sqrt(σ²sX), sqrt(σ²sY)], γ = [γ₁, γ₂], w = [w₁, w₂]),
                       (slope = slope, intercept = intercept))
end
# Nonequivalent Goups : Chained linear Observed Score Equating
struct ResultChainedLinear <: NEATEquateMethod
    table::DataFrame
    synthetic::DataFrame
    estimates::NamedTuple
end
function ChainedLinear(X::NEAT, Y::NEAT)
    # ******************************************** #
    # 1. put X on the scale of V -call lV(x);
    # 2. put V on the scale of Y - call lY(v);
    # 3. obtain Y-equivalent as lY(x).
    # ******************************************** #
    # test score
    x = X.rawX; xv = X.rawV
    y = Y.rawX; yv = Y.rawV
    # ObservableStats
    μx, σx, μxv, σxv, covxv, corxv = ObservableStats(X)
    μy, σy, μyv, σyv, covyv, coryv = ObservableStats(Y)
    # regression slope
    γ₁ = σx / σxv
    γ₂ = σy / σyv
    # estimate
    slope = (σy/σyv)/(σx/σxv)
    intercept = μy + σy/σyv *(μxv - μyv) - slope * μx
    lYx = @. X.tabX.scale * slope + intercept
    tbl =  DataFrame(scaleX = X.tabX.scale, lYx = lYx)
    ResultChainedLinear(tbl,
                        DataFrame(Group = [1,2], γ = [γ₁, γ₂]),
                        (slope = slope, intercept = intercept))
end
# Nonequivalent Goups : Frequency Estimation
struct ResultFrequencyEstimation <: NEATEquateMethod
    table::DataFrame
    marginalX::Matrix
    marginalY::Matrix
end
function FrequencyEstimation(X::NEAT, Y::NEAT; w₁ = length(X.rawX) / (length(X.rawX) + length(Y.rawX)), w₂ = 1.0 - w₁)
    # synthetic weight
    W = w₁ + w₂
    w₁ = w₁ / W; w₂ = w₂ / W
    # prior (the weights from common part)
    J = length(X.tabV.freq)
    h₁ = X.tabV.freq / sum(X.tabV.freq)
    h₂ = Y.tabV.freq / sum(Y.tabV.freq)
    f₂x = zeros(Float64, length(X.tabX.freq))
    g₁y = zeros(Float64, length(Y.tabX.freq))
    for j in 1:length(X.tabX.scale)
        f₂x[j] = X.marginal[j,:]' * h₂
    end
    for j in 1:length(Y.tabX.scale)
        g₁y[j] = Y.marginal[j,:]' * h₁
    end
    # synthetic population
    fsx = @. w₁ * X.tabX.freq + w₂ * f₂x
    fsy = @. w₁ * g₁y + w₂ * Y.tabX.freq
    # Equipercentile Equating
    ftX = FreqTab(DataFrame(scale = X.tabX.scale, freq = fsx, cumprob = cumsum(fsx) ./ sum(fsx)),
                  X.rawX, X.intervalX)
    ftY = FreqTab(DataFrame(scale = Y.tabX.scale, freq = fsy, cumprob = cumsum(fsy) ./ sum(fsy)),
                  Y.rawX, Y.intervalX)
    tbl = equipercentile(ftX, ftY)
    return ResultFrequencyEstimation(tbl, X.marginal, Y.marginal)
end
# Nonequivalent Groups : Braun-Holland Linear Method
struct ResultBraunHolland <: NEATEquateMethod
    table::DataFrame
    marginalX::Matrix
    marginalY::Matrix
end
function BraunHolland(X::NEAT, Y::NEAT; w₁ = length(X.rawX) / (length(X.rawX) + length(Y.rawX)), w₂ = 1.0 - w₁)
    # synthetic weight
    W = w₁ + w₂
    w₁ = w₁ / W; w₂ = w₂ / W
    # prior (the weights from common part)
    J = length(X.tabV.freq)
    h₁ = X.tabV.freq / sum(X.tabV.freq)
    h₂ = Y.tabV.freq / sum(Y.tabV.freq)
    f₂x = zeros(Float64, length(X.tabX.freq))
    g₁y = zeros(Float64, length(Y.tabX.freq))
    for j in 1:length(X.tabX.scale)
        f₂x[j] = X.marginal[j,:]' * h₂
    end
    for j in 1:length(Y.tabX.scale)
        g₁y[j] = Y.marginal[j,:]' * h₁
    end
    # synthetic population
    fsx = @. w₁ * X.tabX.freq + w₂ * f₂x; fsx = fsx ./ sum(fsx)
    fsy = @. w₁ * g₁y + w₂ * Y.tabX.freq; fsy = fsy ./ sum(fsy)
    # synthetic pupulation parameter
    μsx = fsx' * X.tabX.scale; σsx = sqrt(fsx' * (X.tabX.scale .- μsx) .^2)
    μsy = fsy' * Y.tabX.scale; σsy = sqrt(fsy' * (Y.tabX.scale .- μsy) .^2)
    μxv = mean(X.rawV); σxv = std(X.rawV)
    μyv = mean(Y.rawV); σyv = std(Y.rawV)
    # internal regression parameter
    γ₁ = σsx / σxv; γ₂ = σsy / σyv
    # external regression parameter
    slope = γ₂ / γ₁; intercept = μsy + σsy/σyv*(μxv-μyv) - slope * μsx
    lYx = @. X.tabX.scale * slope + intercept
    tbl = DataFrame(scaleX = X.tabX.scale, lYx = lYx)
    return ResultBraunHolland(tbl, X.marginal, Y.marginal)
end

# Nonequivalent Groups : Chained Equipercentile Method
struct ResultChainedEquipercentile <: NEATEquateMethod
    table::DataFrame
end
function ChainedEquipercentile(X::NEAT, Y::NEAT; case = :middle)
    #
    ftX = freqtab(X.rawX); ftXV = freqtab(X.rawV)
    eV₁ = Equipercentile(ftX, ftXV)
    eY₂ = Equipercentile(freqtab(Y.rawV), freqtab(Y.rawX))
    # Search percentile of score V on scale Y
    eYxu = zeros(Float64, length(eY₂.table.scaleY)); eYxl = zeros(Float64, length(eY₂.table.scaleY))
    for (i,v) in enumerate(eY₂.table.eYx)
        P = PRF(v, ftXV)
        eYxu[i] = PFu(P, ftX)
        eYxl[i] = PFl(P, ftX)
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
    tbl = DataFrame(scaleY = eY₂.table.scaleY, eYx = eYx)
    return ResultChainedEquipercentile(tbl)
end
#-----------------
end # module

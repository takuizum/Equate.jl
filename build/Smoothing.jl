# LogLinear Transformation

"""
    LogLinearFormula(df::Int64)

Returns GLM formula can be used in `glm` with `FreqTab` data frame.
`df`, is abbreviation for degrees of formula, represents degree of polynomial log-linear method.
This function works very slowly. If a fixed degree formula will be used repeatedly, define the formula as a variable.
The similer macro `@LogLinearFormula` is provided and it is faster than the function.
"""
function LogLinearFormula(df)
    fml = "@formula(freq ~ 1 +"
    @simd for d in 1:df
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
    fit
end
"""
    presmoothing(F::FreqTab, @LogLinearFormula(df::Int64))

Returns presmoothed frequency table as `SmoothedFreqTab` and `glm` fitted object.

Preserving first C moments of original frequency data, passed `LogLinearFormula(C)` to `fml`.
C is a degree of freedom
"""
function presmoothing(F::EG, fml)
    fit1 = glm(fml, F.tab, Poisson(), LogLink())
    freq = predict(fit1, DataFrame(scale = F.tab.scale))
    tab = DataFrame(scale = F.tab.scale, freq = freq, cumfreq = cumsum(freq),
                    prob = freq ./ sum(freq), cumprob = cumsum(freq) ./ sum(freq))
    return SmoothedFreqTab(tab, F.raw, F.interval, fit1)
end

# function presmoothing(F::EG, degree::Int64)
#     AIC, BIC, AICC, DEVIANCE, LL = zeros.(Float64, fill(degree, 5))
#     for d in 1:degree
#         println("Fitting dof = $degree model.")
#         fml = LogLinearFormula(d)
#         fit = glm(fml, F.tab, Poisson(), LogLink())
#         AIC[d] = aic(fit)
#         BIC[d] = bic(fit)
#         AICC[d] = aicc(fit)
#         LL[d] = loglikelihood(fit)
#         DEVIANCE[d] = deviance(fit)
#     end
#     DataFrame(degree = [1:1:degree;], aic = AIC, aicc = IACC, bic = BIC, loglik = LL, deviance = DEVIANCE)
# end

struct SmoothedNEATFreqTab <: NEAT
    tabX
    tabV
    rawX
    rawV
    intervalX
    intervalV
    fitX
    fitV
    marginal
end
"""
    presmoothing(F::NEAT, LogLinearFormula(df::Int64), LogLinearFormula(df::Int64))

Returns presmoothed frequency table as `SmoothedNEATFreqTab` and `glm` fitted object, 
but `interval` and `marginal` elements, which are returned, are not based on smoothed score.

Preserving first C moments of original frequency data, passed `LogLinearFormula(C)` to `fml`.
C is a degree of freedom
"""
function presmoothing(F::NEAT, fmlX, fmlV)
    # Smoothing X (independent part)
    fitX = glm(fmlX, F.tabX, Poisson(), LogLink())
    freqX = predict(fitX, DataFrame(scale = F.tabX.scale))
    tabX = DataFrame(scale = F.tabX.scale, freq = freqX, cumfreq = cumsum(freqX),
                    prob = freqX ./ sum(freqX), cumprob = cumsum(freqX) ./ sum(freqX))
    # Smoothing V (common part)
    fitV = glm(fmlV, F.tabV, Poisson(), LogLink())
    freqV = predict(fitV, DataFrame(scale = F.tabV.scale))
    tabV = DataFrame(scale = F.tabV.scale, freq = freqV, cumfreq = cumsum(freqV),
                    prob = freqV ./ sum(freqV), cumprob = cumsum(freqV) ./ sum(freqV))
    return SmoothedNEATFreqTab(tabX, tabV, F.rawX, F.rawV, F.intervalX, F.intervalV, fitX, fitV, F.marginal)
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
"""
    KernelSmoothing(X::SmoothedFreqTab; kernel = :Gaussian, hX = 0.622, scale = X.tab.scale)

Returns Kernel smoothed frequency table as `KernelFreqTab`.

# Arguments

- `X` the object which class is `freqtab`

## Optional Arguments

- `kernel`
- `hX` The band width.
- `newint` A new interval to depict the continuized probability line.

# Value

`KernelFreqTab`

In the kernel smoothing, the choice of bandwidth `hX` is an important consideration. When bandwidth is large, the smoothed distribution becomes linear equating fucnction. Contrary, bandwidth is small, it becomes spikye shaped distribution. We recommend to use default value `hX = 6.22` or select it by using `EstBandwidth`.
"""
function KernelSmoothing(X::EG; kernel = :Gaussian, hX = 0.622, newint = X.interval/100)
    scale = (minimum(X.tab.scale) - X.interval/2):newint:(maximum(X.tab.scale) + X.interval/2)
    # hX = bandwidth of cumulative distribution function
    N = sum(X.tab.freq)
    μ = mean(X.raw); σ² = var(X.raw)
    a² = σ² / (σ² + hX^2)
    a = sqrt(a²)
    𝒇hX = zeros(Float64, length(scale))
    FhX = zeros(Float64, length(scale))
    for (i, x) in enumerate(scale), (j, xⱼ) in enumerate(X.tab.scale)
        𝒇hX[i] += X.tab.prob[j]*pdf.(Normal(0, 1), RjX(x, xⱼ, a, μ, hX)) / (a*hX)
        FhX[i] += X.tab.prob[j]*cdf.(Normal(0, 1), RjX(x, xⱼ, a, μ, hX))
    end
    𝒇hX = 𝒇hX ./ sum(𝒇hX) # normalize
    tbl = DataFrame(scale = scale, freq = N .* 𝒇hX, cunfreq = cumsum(N .* 𝒇hX), prob = 𝒇hX, cumprob = cumsum(𝒇hX))
    return KernelFreqTab(tbl, X.raw, X.interval, hX)
end

function BandwidthPenalty(hX, X::SmoothedFreqTab, lprob, rprob, K; kernel = :Gaussian)
    hX = exp(hX[1])
    r = X.tab.prob
    μ = mean(X.raw); σ² = var(X.raw)
    a² = σ² / (σ² + hX^2)
    a =sqrt(a²)
    𝒇hX = zeros(Float64, length(X.tab.scale))
    r𝒇′hX = zeros(Float64, length(X.tab.scale))
    l𝒇′hX = zeros(Float64, length(X.tab.scale))
    for (i, x) in enumerate(X.tab.scale), (j, xⱼ) in enumerate(X.tab.scale)
        R = RjX(x, xⱼ, a, μ, hX)
        𝒇hX[i] += X.tab.prob[j]*pdf.(Normal(0, 1), R) / (a*hX)
        r𝒇′hX[i] -= rprob[j]*pdf.(Normal(0, 1), R) / (a*hX)^2 * R
        l𝒇′hX[i] -= lprob[j]*pdf.(Normal(0, 1), R) / (a*hX)^2 * R
    end
    pen1 = sum((r .- 𝒇hX) .^2)
    pen2 = sum((l𝒇′hX .< 0) .* (r𝒇′hX .> 0))
    return pen1 + K * pen2
end

"""
    EstBandwidth(X::EG; kernel = :Gaussian)

Estimate optimal bandwidth `hX` in Gaussian kernel.

# Arguments

- `X` is the object which the class is `freqtab`.

# Value

Optimize information.

"""
function EstBandwidth(X::SmoothedFreqTab; kernel = :Gaussian, K = 1)
    x = X.tab.scale
    lscore = x .- X.interval/4
    rscore = x .+ X.interval/4
    lpred = predict(X.fit, DataFrame(scale = lscore))
    rpred = predict(X.fit, DataFrame(scale = rscore))
    opt = optimize(hX -> BandwidthPenalty(hX, X, lpred, rpred, K; kernel = kernel), [0.5], method = BFGS())
    println("Minimizer sould be transformed `exp()` before interpretation. \nMinimizer (after taking `exp()`) = $(exp(opt.minimizer[1]))\n\n")
    return opt
end

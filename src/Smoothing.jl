# LogLinear Transformation

"""
    LogLinearFormula(df::Int64)

Returns GLM formula can be used in `glm` with `FreqTab` data frame.

`df`, is abbreviation for degrees of formula, represents degree of polynomial log-linear method.
"""
function LogLinearFormula(df::Int64)
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
"""
    presmoothing(F::FreqTab; fml = LogLinearFormula(df::Int64))

Returns presmoothed frequency table as `SmoothedFreqTab` and `glm` fitted object.

Preserving first C moments of original frequency data, passed `LogLinearFormula(C)` to `fml`.
"""
function presmoothing(F::EG, fml)
    fit1 = glm(fml, F.tab, Poisson(), LogLink())
    freq = predict(fit1, DataFrame(scale = F.tab.scale))
    tab = DataFrame(scale = F.tab.scale, freq = freq, cumfreq = cumsum(freq),
                    pred = freq ./ sum(freq), cumprob = cumsum(freq) ./ sum(freq))
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
"""
    KernelSmoothing(X::EG; kernel = :Gaussian, hX = 0.622, scale = X.tab.scale)

Returns Kernel smoothed frequency table as `KernelFreqTab`.

In the kernel smoothing, the choice of bandwidth `hX` is an important consideration. When bandwidth is large, the smoothed distribution becomes linear equating fucnction. Contrary, bandwidth is small, it becomes spikye shaped distribution. We recommend to use default value `hX = 6.22` or select it by using `EstBandwidth`
"""
function KernelSmoothing(X::EG; kernel = :Gaussian, hX = 0.622, scale = X.tab.scale)
    # hX = bandwidht of cumulative distribution function
    μ = mean(X.raw); σ² = var(X.raw)
    a² = σ² / (σ² + hX^2)
    a = sqrt(a²)
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

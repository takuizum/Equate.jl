# LogLinear Transformation

"""
    LogLinearFormula(df::Int64)

Returns GLM formula can be used in `glm` with `FreqTab` data frame.
`df`, is abbreviation for degrees of formula, represents degree of polynomial log-linear method.
This function works very slowly. If a fixed degree formula will be used repeatedly, define the formula as a variable.
"""
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

const fml₁ = @formula(freq ~ 1 + scale)
const fml₂ = @formula(freq ~ 1 + scale + scale^2)
const fml₃ = @formula(freq ~ 1 + scale + scale^2 + scale^3)
const fml₄ = @formula(freq ~ 1 + scale + scale^2 + scale^3 + scale^4)
const fml₅ = @formula(freq ~ 1 + scale + scale^2 + scale^3 + scale^4 + scale^5)
const fml₆ = @formula(freq ~ 1 + scale + scale^2 + scale^3 + scale^4 + scale^5 + scale^6)
const fml₇ = @formula(freq ~ 1 + scale + scale^2 + scale^3 + scale^4 + scale^5 + scale^6 + scale^7)
const fml₈ = @formula(freq ~ 1 + scale + scale^2 + scale^3 + scale^4 + scale^5 + scale^6 + scale^7 + scale^8)
const fml₉ = @formula(freq ~ 1 + scale + scale^2 + scale^3 + scale^4 + scale^5 + scale^6 + scale^7 + scale^8 + scale^9)
const fml₁₀ = @formula(freq ~ 1 + scale + scale^2 + scale^3 + scale^4 + scale^5 + scale^6 + scale^7 + scale^8 + scale^9 + scale^10)

mutable struct SmoothedSGFreqTab <: EG
    table
    raw
    interval
    fit
    method
    model
end

mutable struct SmoothedNEATFreqTab <: NEAT
    tableX
    tableV
    rawX # independent form
    rawV # common form
    intervalX
    intervalV
    marginal # conditional freqency
    statsX
    statsV
    method
    model
end
"""
    presmoothing(F::EG, fml)
    presmoothing(F::EG, UpToDgree::Int64)
    presmoothing(F::NEAT, fml)


Returns presmoothed frequency table as `SmoothedSG(NEAT)FreqTab` and `glm` fitted object.

Preserving first C moments of original frequency data. `LogLinearFormula(C)` creates `fml`.

C is a degree of freedom.

For more instance, use, `Equate.fml₆` to fit the log linear model with 6 dof.

# Examples
```julia
julia> tab = freqtab(expandtable(ACTmath.scale, ACTmath.xcount));
julia> presmoothing(tab, Equate.fml₄)
Equate.SmoothedFreqTab(40×5 DataFrame
│ Row │ scale   │ freq     │ cumfreq  │ prob        │ cumprob     │
│     │ Float64 │ Float64? │ Float64? │ Float64     │ Float64     │
├─────┼─────────┼──────────┼──────────┼─────────────┼─────────────┤
│ 1   │ 1.0     │ 0.864929 │ 0.864929 │ 0.000199799 │ 0.000199799 │
│ 2   │ 2.0     │ 2.47529  │ 3.34022  │ 0.000571794 │ 0.000771592 │
⋮
│ 38  │ 38.0    │ 31.9351  │ 4291.07  │ 0.00737703  │ 0.991238    │
│ 39  │ 39.0    │ 22.7842  │ 4313.85  │ 0.00526315  │ 0.996501    │
│ 40  │ 40.0    │ 15.1484  │ 4329.0   │ 0.00349928  │ 1.0         │, [1, 2, 3, 3, 3, 4, 4, 4, 4, 4  …  40, 40, 40, 40, 40, 40, 40, 40, 40, 40], 1.0, StatsModels.TableRegressionModel{GLM.GeneralizedLinearModel{GLM.GlmResp{Array{Float64,1},Distributions.Poisson{Float64},GLM.LogLink},GLM.DensePredChol{Float64,LinearAlgebra.Cholesky{Float64,Array{Float64,2}}}},Array{Float64,2}}

freq ~ 1 + scale + :(scale ^ 2) + :(scale ^ 3) + :(scale ^ 4)

Coefficients:
─────────────────────────────────────────────────────────────────────────────────
                   Coef.   Std. Error       z  Pr(>|z|)    Lower 95%    Upper 95%
─────────────────────────────────────────────────────────────────────────────────
(Intercept)  -1.35982     0.310243      -4.38    <1e-4   -1.96789     -0.751757
scale         1.30127     0.0728808     17.85    <1e-70   1.15843      1.44412
scale ^ 2    -0.089081    0.00589528   -15.11    <1e-50  -0.100636    -0.0775265
scale ^ 3     0.00254824  0.000195003   13.07    <1e-38   0.00216604   0.00293044
scale ^ 4    -2.67698e-5  2.25114e-6   -11.89    <1e-31  -3.1182e-5   -2.23576e-5
─────────────────────────────────────────────────────────────────────────────────)

```
"""
function presmoothing(F::EG, fml)
    fit1 = glm(fml, F.table, Poisson(), LogLink())
    freq = predict(fit1, DataFrame(scale = F.table.scale))
    tab = DataFrame(scale = F.table.scale, freq = freq, cumfreq = cumsum(freq),
                    prob = freq ./ sum(freq), cumprob = cumsum(freq) ./ sum(freq))
    return SmoothedSGFreqTab(tab, F.raw, F.interval, fit1, :loglinear, fml)
end

function presmoothing(F::NEAT, fml)
    longtab = DataFrame(
        s = reshape(F.marginal, (*(size(F.marginal)...), 1))[:],
        x = vcat(fill(F.tableX.scale, length(F.tableV.scale))...),
        v = vcat(fill.(F.tableV.scale, length(F.tableX.scale))...)
    )
    fit1 = glm(@formula(s ~ x^1*v^1 + x^2*v^2), longtab, Poisson(), LogLink())
    freq = predict(fit1, longtab)
    marginalfreq = reshape(freq, size(F.marginal))
    freqX = sum(marginalfreq, dims = 2)
    freqV = sum(marginalfreq, dims = 1)[:]
    tabX = DataFrame(scale = F.tableX.scale, freq = freqX, cumfreq = cumsum(freqX),
                     prob = freqX ./ sum(freqX), cumprob = cumsum(freqX) ./ sum(freqX))
    tabV = DataFrame(scale = F.tableV.scale, freq = freqV, cumfreq = cumsum(freqV),
                     prob = freqV ./ sum(freqV), cumprob = cumsum(freqV) ./ sum(freqV))
    return SmoothedNEATFreqTab(
        tabX, tabV,
        F.rawX, F.rawV,
        F.intervalX, F.intervalV,
        marginalfreq, 
        basicstats(tabX, F.rawX), basicstats(tabV, F.rawV)
    )
end

function basicstats(F::DataFrame, X)
    N = sum(F.freq)
    Nm = sum(X .=== missing)
    mins, maxs = minimum(X), maximum(X)
    μ = F.scale'F.freq
    σ = sqrt(sum((F.scale .- μ).^2 .* F.freq))
    k = sum((F.scale .- μ).^3 .* F.prob) / σ^3
    s = sum((F.scale .- μ).^4 .* F.prob) / σ^4 - 3
    (N = N, Missing = Nm, min = mins, max = maxs, μ = μ, σ = σ, kurtosis = k, skewness = s)
end

"""
    presmoothing(F::EG, degree::Int64)
Examine various degree of freedom for fitting the log liner models to the raw score histogram.

# Examples
```julia
julia> tab = freqtab(expandtable(ACTmath.scale, ACTmath.xcount))
Equate.FreqTab(40×5 DataFrame
│ Row │ scale   │ freq  │ cumfreq │ prob       │ cumprob  │
│     │ Float64 │ Int64 │ Int64   │ Float64    │ Float64  │
├─────┼─────────┼───────┼─────────┼────────────┼──────────┤
│ 1   │ 1.0     │ 1     │ 1       │ 0.000231   │ 0.000231 │
│ 2   │ 2.0     │ 1     │ 2       │ 0.000231   │ 0.000462 │
⋮
│ 38  │ 38.0    │ 38    │ 4291    │ 0.00877801 │ 0.991222 │
│ 39  │ 39.0    │ 23    │ 4314    │ 0.00531301 │ 0.996535 │
│ 40  │ 40.0    │ 15    │ 4329    │ 0.003465   │ 1.0      │, [1, 2, 3, 3, 3, 4, 4, 4, 4, 4  …  40, 40, 40, 40, 40, 40, 40, 40, 40, 40], 1.0)

julia> fit = presmoothing(tab, 10)
Fitting dof = 1 model.
Fitting dof = 2 model.
Fitting dof = 3 model.
Fitting dof = 4 model.
Fitting dof = 5 model.
Fitting dof = 6 model.
Fitting dof = 7 model.
Fitting dof = 8 model.
Fitting dof = 9 model.
Fitting dof = 10 model.
10×7 DataFrame. Omitted printing of 1 columns
│ Row │ degree │ aic     │ aicc    │ bic     │ loglik   │ deviance │
│     │ Int64  │ Float64 │ Float64 │ Float64 │ Float64  │ Float64  │
├─────┼────────┼─────────┼─────────┼─────────┼──────────┼──────────┤
│ 1   │ 1      │ 2235.18 │ 2235.5  │ 2238.56 │ -1115.59 │ 1988.28  │
│ 2   │ 2      │ 661.558 │ 662.225 │ 666.625 │ -327.779 │ 412.654  │
⋮
│ 8   │ 8      │ 290.675 │ 296.675 │ 305.875 │ -136.338 │ 29.7704  │
│ 9   │ 9      │ 292.573 │ 300.159 │ 309.462 │ -136.286 │ 29.6683  │
│ 10  │ 10     │ 294.051 │ 303.48  │ 312.629 │ -136.026 │ 29.1466  │

```
"""
function presmoothing(F::EG, degree::Int64)
    if degree > 11
        throw(ArgumentError("Maximum degree for examining must be less than ore equal to 10."))
    end
    fmls = (fml₁, fml₂, fml₃, fml₄, fml₅, fml₆, fml₇, fml₈, fml₉, fml₁₀)
    AIC, BIC, AICC, DEVIANCE, LL = zeros.(Float64, fill(degree, 5))
    FIT = Any[]
    for d in 1:degree
        println("Fitting dof = $d model.")
        fit = glm(fmls[d], F.table, Poisson(), LogLink())
        push!(FIT, fit)
        AIC[d] = StatsBase.aic(fit)
        BIC[d] = StatsBase.bic(fit)
        AICC[d] = StatsBase.aicc(fit)
        LL[d] = StatsBase.loglikelihood(fit)
        DEVIANCE[d] = StatsBase.deviance(fit)
    end
    DataFrame(degree = [1:1:degree;], aic = AIC, aicc = AICC, bic = BIC, loglik = LL, deviance = DEVIANCE, fit = FIT)
end

"""
    presmoothing(F::NEAT, LogLinearFormula(df::Int64), LogLinearFormula(df::Int64))

Returns presmoothed frequency table as `SmoothedNEATFreqTab`, but `interval` and `marginal` elements are not based on smoothed score.
This method fit the loglinear model for individual, unidimensional frequency table.

Preserving first C moments of original frequency data, passed `LogLinearFormula(C)` to `fml`.
C is a degree of freedom
"""
function presmoothing(F::NEAT, fmlX, fmlV)
    # Smoothing X (independent part)
    fitX = glm(fmlX, F.tableX, Poisson(), LogLink())
    freqX = predict(fitX, DataFrame(scale = F.tableX.scale))
    tabX = DataFrame(scale = F.tableX.scale, freq = freqX, cumfreq = cumsum(freqX),
                    prob = freqX ./ sum(freqX), cumprob = cumsum(freqX) ./ sum(freqX))
    # Smoothing V (common part)
    fitV = glm(fmlV, F.tableV, Poisson(), LogLink())
    freqV = predict(fitV, DataFrame(scale = F.tableV.scale))
    tabV = DataFrame(scale = F.tableV.scale, freq = freqV, cumfreq = cumsum(freqV),
                    prob = freqV ./ sum(freqV), cumprob = cumsum(freqV) ./ sum(freqV))
    return SmoothedNEATFreqTab(tabX, tabV, F.rawX, F.rawV, F.intervalX, F.intervalV, fitX, fitV, F.marginal)
end


# Kernel method
function RjX(x, xⱼ, a, μ, hX)
    return (x - a * xⱼ - (1-a)*μ) / (a*hX)
end
struct KernelFreqTab <: EG
    table
    raw
    interval
    Bandwidth
end
"""
    KernelSmoothing(X::EG; kernel = :Gaussian, hX = 0.622, scale = X.table.scale)

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
    scale = (minimum(X.table.scale) - X.interval/2):newint:(maximum(X.table.scale) + X.interval/2)
    # hX = bandwidth of cumulative distribution function
    N = sum(X.table.freq)
    μ = mean(X.raw); σ² = var(X.raw)
    a² = σ² / (σ² + hX^2)
    a = sqrt(a²)
    𝒇hX = zeros(Float64, length(scale))
    FhX = zeros(Float64, length(scale))
    for (i, x) in enumerate(scale), (j, xⱼ) in enumerate(X.table.scale)
        𝒇hX[i] += X.table.prob[j]*pdf.(Normal(0, 1), RjX(x, xⱼ, a, μ, hX)) / (a*hX)
        FhX[i] += X.table.prob[j]*cdf.(Normal(0, 1), RjX(x, xⱼ, a, μ, hX))
    end
    𝒇hX = 𝒇hX ./ sum(𝒇hX) # normalize
    tbl = DataFrame(scale = scale, freq = N .* 𝒇hX, cunfreq = cumsum(N .* 𝒇hX), prob = 𝒇hX, cumprob = cumsum(𝒇hX))
    return KernelFreqTab(tbl, X.raw, X.interval, hX)
end

function BandwidthPenalty(hX, X::EG, lprob, rprob, K; kernel = :Gaussian)
    hX = exp(hX[1])
    r = X.table.prob
    μ = mean(X.raw); σ² = var(X.raw)
    a² = σ² / (σ² + hX^2)
    a =sqrt(a²)
    𝒇hX = zeros(Float64, length(X.table.scale))
    r𝒇′hX = zeros(Float64, length(X.table.scale))
    l𝒇′hX = zeros(Float64, length(X.table.scale))
    for (i, x) in enumerate(X.table.scale), (j, xⱼ) in enumerate(X.table.scale)
        R = RjX(x, xⱼ, a, μ, hX)
        𝒇hX[i] += X.table.prob[j]*pdf.(Normal(0, 1), R) / (a*hX)
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
function EstBandwidth(X::EG; kernel = :Gaussian, K = 1)
    x = X.table.scale
    lscore = x .- X.interval/4
    rscore = x .+ X.interval/4
    lpred = predict(X.fit, DataFrame(scale = lscore))
    rpred = predict(X.fit, DataFrame(scale = rscore))
    opt = optimize(hX -> BandwidthPenalty(hX, X, lpred, rpred, K; kernel = kernel), [0.5], method = BFGS())
    println("Minimizer sould be transformed `exp()` before interpretation. \nMinimizer (after taking `exp()`) = $(exp(opt.minimizer[1]))\n\n")
    return opt
end

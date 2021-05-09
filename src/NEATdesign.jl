# Nonequivalent Groups : Linear methods
function ObservableStats(F::NEAT)
    x = F.rawX; v = F.rawV
    μx = mean(x); σx = std(x)
    μv = mean(v); σv = std(v)
    covxv = cov(x, v); corxv = cor(x, v)
    return μx, σx, μv, σv, covxv, corxv
end

function SummaryStats(F::NEAT)
    r(x) = round2(x; digits = 2)
    x = F.rawX; v = F.rawV
    μx = mean(x); σx = std(x); kx = kurtosis(x); sx = skewness(x); 
    μv = mean(v); σv = std(v); kv = kurtosis(v); sv = skewness(v);
    df = DataFrame(
        test = ["X", "V"], mean = [μx, μv], sigma = [σx, σv], kurtosis = [kx, kv], skewness = [sx, sv],
        min = [minimum(x), minimum(v)], max = [maximum(x), maximum(v)], N = [length(x), length(v)]
    )
    println(df)
    return df;
end

mutable struct NEATEquateResult <: NEATEquateMethod
    method
    table
    synthetic
    estimates
    data
end

"""
    Tucker(X::NEAT, Y::NEAT; w₁ = length(X.rawX) / (length(X.rawX) + length(Y.rawX)), w₂ = 1.0 - w₁)

Tucker observed score method.
"""
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
    lYx = @. X.tableX.scale * slope + intercept
    tbl = DataFrame(scaleX = X.tableX.scale, lYx = lYx)
    return NEATEquateResult(
        :Tucker,
        tbl,
        DataFrame(Group = [1, 2], μ = [μsX, μsY], σ = [sqrt(σ²sX), sqrt(σ²sY)], γ = [γ₁, γ₂], w = [w₁, w₂]),
        (slope = slope, intercept = intercept), 
        (X = X, Y = Y)
    )
end
# Nonequivalent Groups : Levine under a classical congeneric model
"""
    LevineCongeneric(X::NEAT, Y::NEAT; common = :external, w₁ = length(X.rawX) / (length(X.rawX) + length(Y.rawX)), w₂ = 1.0 - w₁)

Levine methods under a classical congeneric model.
"""
function LevineCongeneric(X::NEAT, Y::NEAT; common = :external, w₁ = length(X.rawX) / (length(X.rawX) + length(Y.rawX)), w₂ = 1.0 - w₁)
    W = w₁ + w₂
    w₁ = w₁ / W; w₂ = w₂ / W
    # test score
    x = X.rawX; xv = X.rawV
    y = Y.rawX; yv = Y.rawV
    # ObservableStats
    μx, σx, μxv, σxv, covxv, corxv = ObservableStats(X)
    μy, σy, μyv, σyv, covyv, coryv = ObservableStats(Y)
    # regression slope
    if common == :internal
        γ₁ = σx^2 / covxv
        γ₂ = σy^2 / covyv
    elseif common == :external
        γ₁ = (σx^2 + covxv) / (σxv^2 + covxv)
        γ₂ = (σy^2 + covyv) / (σyv^2 + covyv)
    else
        println("`common` argument can recieve only :internal or :external")
    end
    # synthetic mean and var
    μsX = μx - w₂*γ₁*(μxv-μyv)
    μsY = μy + w₁*γ₂*(μxv-μyv)
    σ²sX = σx^2 - w₂*γ₁^2*(σxv^2-σyv^2) + w₁*w₂*γ₁^2*(μxv-μyv)^2
    σ²sY = σy^2 + w₁*γ₂^2*(σxv^2-σyv^2) + w₁*w₂*γ₂^2*(μxv-μyv)^2
    # transformation
    slope = sqrt(σ²sY)/sqrt(σ²sX); intercept = μsY - slope*μsX
    lYx = @. X.tableX.scale * slope + intercept
    tbl = DataFrame(scaleX = X.tableX.scale, lYx = lYx)
    return NEATEquateResult(
        :LevineCongeneric, 
        tbl,
        DataFrame(Group = [1, 2], μ = [μsX, μsY], σ = [sqrt(σ²sX), sqrt(σ²sY)], γ = [γ₁, γ₂], w = [w₁, w₂]),
        (slope = slope, intercept = intercept), 
        (X = X, Y = Y)
    )
end
# Nonequivalent Groups : Chained linear Observed Score Equating
"""
    ChainedLinear(X::NEAT, Y::NEAT)

Equate 2 tests in a chain through the common test scores. Now supports only 2 form equating.

# Basic Algorithm

1. put X on the scale of V -call lV(x);
2. put V on the scale of Y - call lY(v);
3. obtain Y-equivalent as lY(x).


"""
function ChainedLinear(X::NEAT, Y::NEAT)
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
    slope = γ₂/γ₁
    intercept = μy + γ₂ *(μxv - μyv) - slope * μx
    lYx = @. X.tableX.scale * slope + intercept
    tbl =  DataFrame(scaleX = X.tableX.scale, lYx = lYx)
    return NEATEquateResult(
        :ChainedLinear, 
        tbl,
        DataFrame(Group = [1,2], γ = [γ₁, γ₂]),
        (slope = slope, intercept = intercept), 
        (X = X, Y = Y)
    )
end

"""
    ChainedMean(X::NEAT, Y::NEAT)
Chained mean equating, which results are identical for Levine's observed score and true score (not implemented) when w₁ = 1.0 (there is no option for the population weights in the mean equating).
"""
function ChainedMean(X::NEAT, Y::NEAT)
    # test score
    x = X.rawX; xv = X.rawV
    y = Y.rawX; yv = Y.rawV
    # ObservableStats
    μx, σx, μxv, σxv, covxv, corxv = ObservableStats(X)
    μy, σy, μyv, σyv, covyv, coryv = ObservableStats(Y)
    γ₂ = σy / σyv
    # estimate
    slope = 1.0
    intercept = μy + γ₂ *(μxv - μyv) - slope * μx
    lYx = @. X.tableX.scale * slope + intercept
    tbl =  DataFrame(scaleX = X.tableX.scale, lYx = lYx)
    return NEATEquateResult(
        :ChainedMean, 
        tbl,
        DataFrame(Group = [1,2], γ = [γ₁, γ₂]),
        (slope = slope, intercept = intercept), 
        (X = X, Y = Y)
    )
end
# Nonequivalent Groups : Braun-Holland Linear Method
"""
    BraunHolland(X::NEAT, Y::NEAT; w₁ = length(X.rawX) / (length(X.rawX) + length(Y.rawX)), w₂ = 1.0 - w₁)

Conduct Braun-Holland Linear equating, which use means and standard deviations estimated by the synsetic population distributions.

# Arguments

- `X` frequency table generated by `freqtab(X, V)`, `V` is the test score of common item part
- `Y` same as above.
- `w₁`, `w₂` The synsetic population weights, that are constrained to be w₁ + w₂ = 1.
"""
function BraunHolland(X::NEAT, Y::NEAT; w₁ = length(X.rawX) / (length(X.rawX) + length(Y.rawX)), w₂ = 1.0 - w₁)
    # synthetic weight
    W = w₁ + w₂
    w₁ = w₁ / W; w₂ = w₂ / W
    # prior (the weights from common part)
    J = length(X.tableV.freq)
    h₁ = X.tableV.prob
    h₂ = Y.tableV.prob
    # Conditional distribution
    conditionalX = X.marginal ./ sum(X.marginal)
    conditionalY = Y.marginal ./ sum(Y.marginal)
    for v in 1:length(X.tableV.scale)
        conditionalX[:, v] ./= h₁[v]
    end
    for v in 1:length(Y.tableV.scale)
        conditionalY[:, v] ./= h₂[v]
    end
    # Marginal out
    f₂x = zeros(Float64, length(X.tableX.freq))
    g₁y = zeros(Float64, length(Y.tableX.freq))
    for j in 1:length(X.tableX.scale)
        f₂x[j] = conditionalX[j,:]'h₂
    end
    for j in 1:length(Y.tableX.scale)
        g₁y[j] = conditionalY[j,:]'h₁
    end
    # normalize
    f₂x = f₂x / sum(f₂x)
    g₁y = g₁y / sum(g₁y)
    # synthetic population
    fsx = @. w₁ * X.tableX.prob + w₂ * f₂x; fsx = fsx ./ sum(fsx)
    fsy = @. w₁ * g₁y + w₂ * Y.tableX.prob; fsy = fsy ./ sum(fsy)
    @show f₂x, g₁y
    # synthetic pupulation parameter
    μsx = fsx' * X.tableX.scale
    σsx = sqrt(fsx' * (X.tableX.scale .- μsx) .^2)
    μsy = fsy' * Y.tableX.scale
    σsy = sqrt(fsy' * (Y.tableX.scale .- μsy) .^2)
    # external regression parameter
    slope = σsy / σsx; intercept = μsy - slope * μsx
    lYx = @. X.tableX.scale * slope + intercept
    tbl = DataFrame(scaleX = X.tableX.scale, lYx = lYx)
    return NEATEquateResult(
        :BraunHolland, 
        tbl,
        DataFrame(Group = [1, 2], μ = [μsx, μsy], σ = [σsx, σsy], w = [w₁, w₂]),
        (slope = slope, intercept = intercept), 
        (X = X, Y = Y)
    )
end

# Nonequivalent Goups : Frequency Estimation
"""
    FrequencyEstimation(X::NEAT, Y::NEAT; w₁ = length(X.rawX) / (length(X.rawX) + length(Y.rawX)), w₂ = 1.0 - w₁)

Conduct frequency estimation equipercentile equating, which assumes, for both form, the conditional distribution of total score given each common part score is the same in both populations.

"""
function FrequencyEstimation(X::NEAT, Y::NEAT; w₁ = length(X.rawX) / (length(X.rawX) + length(Y.rawX)), w₂ = 1.0 - w₁, case = :middle)
    # synthetic weight
    W = w₁ + w₂
    w₁ = w₁ / W; w₂ = w₂ / W
    # prior (the weights from common part)
    J = length(X.tableV.freq)
    h₁ = X.tableV.prob
    h₂ = Y.tableV.prob
    # Conditional distribution
    conditionalX = X.marginal ./ sum(X.marginal)
    conditionalY = Y.marginal ./ sum(Y.marginal)
    for v in 1:length(X.tableV.scale)
        conditionalX[:, v] ./= h₁[v]
    end
    for v in 1:length(Y.tableV.scale)
        conditionalY[:, v] ./= h₂[v]
    end
    # Marginal out
    f₂x = zeros(Float64, length(X.tableX.freq))
    g₁y = zeros(Float64, length(Y.tableX.freq))
    for j in 1:length(X.tableX.scale)
        f₂x[j] = conditionalX[j,:]'h₂
    end
    for j in 1:length(Y.tableX.scale)
        g₁y[j] = conditionalY[j,:]'h₁
    end
    # normalize
    f₂x = f₂x / sum(f₂x)
    g₁y = g₁y / sum(g₁y)
    # synthetic population
    fsx = @. w₁ * X.tableX.freq + w₂ * f₂x
    fsy = @. w₁ * g₁y + w₂ * Y.tableX.freq
    # Equipercentile Equating
    ftX = FreqTab(DataFrame(scale = X.tableX.scale, freq = fsx, cumprob = cumsum(fsx) ./ sum(fsx)),
                  X.rawX, X.intervalX, (; dummy = "dummy"))
    ftY = FreqTab(DataFrame(scale = Y.tableX.scale, freq = fsy, cumprob = cumsum(fsy) ./ sum(fsy)),
                  Y.rawX, Y.intervalX, (; dummy = "dummy"))
    tbl = Equipercentile(ftX, ftY; case = case)
    return NEATEquateResult(
        :FrequencyEstimation, 
        tbl.table, 
        DataFrame(Group = [1, 2], μ = [mean(fsx), mean(fsy)], σ = [std(fsx), std(fsy)], w = [w₁, w₂]),
        nothing, 
        (X = X, Y = Y)
    )
end

# Nonequivalent Groups : Chained Equipercentile Method
"""
    ChainedEquipercentile(X::NEAT, Y::NEAT; case = :middle)

# Algorithm
Chained equipercentile equating is, in short, to conduct equipercentile equating consecutively.
First, find equipercentile relationship for score X to scores V on the first population.
This function is referred to as 𝑒V1(𝑥).
Second, fubd the equipercentile relationship for converting scores on the common items(V) to scores on the form Y based on examinees from the population 2.
Refer to the resulting function as 𝑒Y2(𝑣).
"""
function ChainedEquipercentile(X::NEAT, Y::NEAT; case = :lower)
    # Find equipercentile relationship for score V on the first population
    ftX = freqtab(X.rawX; interval = X.intervalX, scale = X.tableX.scale)
    ftXV = freqtab(X.rawV; interval = X.intervalV, scale = X.tableV.scale)
    ftYV = freqtab(Y.rawV; interval = Y.intervalV, scale = Y.tableV.scale)
    ftY = freqtab(Y.rawX; interval = Y.intervalX, scale = Y.tableX.scale)
    eV₁x = Equipercentile(ftX, ftXV; case = case)
    # Search percentile of score V on scale Y
    eYxu = zeros(Float64, length(eV₁x.table.scaleX))
    eYxl = zeros(Float64, length(eV₁x.table.scaleX))
    for (i,v) in enumerate(eV₁x.table.eYx)
        P = PRF(v, ftYV)
        eYxu[i] = PFu(P, ftY)
        eYxl[i] = PFl(P, ftY)
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
    tbl = DataFrame(scaleX = eV₁x.table.scaleX, eYx = eYx)
    return NEATEquateResult(
        :ChainedEquipercentile, 
        tbl, 
        nothing, 
        nothing, 
        (X = X, Y = Y)
    )
end
#-----------------

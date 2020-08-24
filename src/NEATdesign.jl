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
    lYx = @. X.tabX.scale * slope + intercept
    tbl = DataFrame(scaleX = X.tabX.scale, lYx = lYx)
    return ResultTucker(tbl,
                       DataFrame(Group = [1, 2], μ = [μsX, μsY], σ = [sqrt(σ²sX), sqrt(σ²sY)], γ = [γ₁, γ₂], w = [w₁, w₂]),
                       (slope = slope, intercept = intercept));
end
# Nonequivalent Groups : Levine under a classical congeneric model
struct ResultLevineCongeneric <: NEATEquateMethod
    table::DataFrame
    synthetic::DataFrame
    estimates::NamedTuple
end
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
    lYx = @. X.tabX.scale * slope + intercept
    tbl = DataFrame(scaleX = X.tabX.scale, lYx = lYx)
    return ResultLevineCongeneric(tbl,
                                  DataFrame(Group = [1, 2], μ = [μsX, μsY], σ = [sqrt(σ²sX), sqrt(σ²sY)], γ = [γ₁, γ₂], w = [w₁, w₂]),
                                  (slope = slope, intercept = intercept))
end
# Nonequivalent Groups : Chained linear Observed Score Equating
struct ResultChainedLinear <: NEATEquateMethod
    table::DataFrame
    synthetic::DataFrame
    estimates::NamedTuple
end
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
"""
    FrequencyEstimation(X::NEAT, Y::NEAT; w₁ = length(X.rawX) / (length(X.rawX) + length(Y.rawX)), w₂ = 1.0 - w₁)

Conduct frequency estimation equipercentile equating, which assumes, for both form, the conditional distribution of total score given each common part score is the same in both populations.

"""
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
    tbl = Equipercentile(ftX, ftY)
    return ResultFrequencyEstimation(tbl.table, X.marginal, Y.marginal)
end
# Nonequivalent Groups : Braun-Holland Linear Method
struct ResultBraunHolland <: NEATEquateMethod
    table::DataFrame
    synthetic::DataFrame
    estimates::NamedTuple
end
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
    return ResultBraunHolland(tbl, 
                              DataFrame(Group = [1, 2], μ = [μsx, μsy], σ = [σsx, σsy], γ = [γ₁, γ₂], w = [w₁, w₂]),
                              (slope = slope, intercept = intercept))
end

# Nonequivalent Groups : Chained Equipercentile Method
struct ResultChainedEquipercentile <: NEATEquateMethod
    table::DataFrame
end
"""
    ChainedEquipercentile(X::NEAT, Y::NEAT; case = :middle)

# Algorithm
Chained equipercentile equating is, in short, to conduct equipercentile equating consecutively.
First, find equipercentile relationship for score X to scores V on the first population.
This function is referred to as 𝑒V1(𝑥).
Second, fubd the equipercentile relationship for converting scores on the common items(V) to scores on the form Y based on examinees from the population 2.
Refer to the resulting function as 𝑒Y2(𝑣).
"""
function ChainedEquipercentile(X::NEAT, Y::NEAT; case = :middle)
    # Find equipercentile relationship for score V on the first population
    ftX = freqtab(X.rawX; interval = X.intervalX, scale = X.tabX.scale)
    ftXV = freqtab(X.rawV; interval = X.intervalV, scale = X.tabV.scale)
    eV₁ = Equipercentile(ftX, ftXV)
    eY₂ = Equipercentile(freqtab(Y.rawV; interval = Y.intervalV, scale = Y.tabV.scale), freqtab(Y.rawX; interval = Y.intervalX, scale = Y.tabX.scale))
    # Search percentile of score V on scale Y
    eYxu = zeros(Float64, length(eV₁.table.scaleX))
    eYxl = zeros(Float64, length(eV₁.table.scaleX))
    for (i,v) in enumerate(eV₁.table.eYx)
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
    tbl = DataFrame(scaleX = eV₁.table.scaleX, eYx = eYx)
    return ResultChainedEquipercentile(tbl)
end
#-----------------

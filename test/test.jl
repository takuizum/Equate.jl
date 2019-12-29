using Distributions, Random

Random.seed!(1234)
X = rand(Poisson(100), 100)
Y = rand(Poisson(90), 100)

# Freqtab
ftX = freqtab(X)
# CDF
CDF(100, ftX)
CDF(0, ftX)

# PRF
t = PRF(120, ftX)
PRF(-1, ftX)
PRF(0, ftX)
PRF(200, ftX)
PRF(101, ftX)

# PF
PFu(100, ftX)
PFu(101, ftX)
PFu(40, ftX)
PFl(40, ftX)
PFl(80, ftX)


# Equipercentile Equating
ftX = freqtab(X); ftY = freqtab(Y)
equipercentile(ftX, ftY; case = :upper) |> print
equipercentile(ftX, ftY; case = :lower) |> print
equipercentile(ftX, ftY; case = :both) |> print
equipercentile(ftX, ftY; case = :middle) |> print

# Log Linear Transformation
fit1 = glm(LogLinearFormula(4), ftX.tab, Poisson(), LogLink())
ftXsmoothed, fit1 = presmoothing(ftX; fml = LogLinearFormula(3))
predict(fit1)
using Plots
plot(ftXsmoothed.tab.scale, predict(fit1), label = "Smoothed degree = 4")
histogram!(X; bins = length(X), alpha = 0.5, label = "observed score")
# link[https://github.com/JuliaStats/StatsModels.jl/blob/master/docs/src/formula.md]

# Linear Equating
linear(ftX, ftY)

# Nonequivalent group design
ftXneg = freqtab(X,Y)
sum(ftXneg.marginal)
ftXneg.marginal |> heatmap

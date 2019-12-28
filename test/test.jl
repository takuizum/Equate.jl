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
PFᵤ(100, ftX)
PFᵤ(101, ftX)
PFᵤ(40, ftX)
PFₗ(40, ftX)
PFₗ(80, ftX)


# Equipercentile Equating
ftX = freqtab(X); ftY = freqtab(Y)
e𝒀x(ftX, ftY; case = :upper) |> print
e𝒀x(ftX, ftY; case = :lower) |> print
e𝒀x(ftX, ftY; case = :both) |> print
e𝒀x(ftX, ftY; case = :middle) |> print

ftX.tab.freq
ftX.tab.scale

fit1 = glm(eval(Meta.parse("@formula(freq ~ scale^1 + scale^2)")), ftX.tab, Poisson(), LogLink())
fit1 = presmoothing(ftX)
predict(fit1)
using Plots
plot(ftX.tab.scale, predict(fit1))
histogram!(X; bins = length(X), alpha = 0.5)
# link[https://github.com/JuliaStats/StatsModels.jl/blob/master/docs/src/formula.md]

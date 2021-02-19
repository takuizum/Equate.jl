include("testdata.jl")
# Equipercentile Equating
using DataFrames
X = fill.(ACTmath.scale, ACTmath.xcount) |> Iterators.flatten |> collect # is euqivalent ExpandTable(ACTmath.scale, ACTmath.xcount)
Y = fill.(ACTmath.scale, ACTmath.ycount) |> Iterators.flatten |> collect
ftX = freqtab(X)
ftY = freqtab(Y)
Equipercentile(ftX, ftY; case = :upper) |> print
Equipercentile(ftX, ftY; case = :lower) |> print
Equipercentile(ftX, ftY; case = :both) |> print
resEp = Equipercentile(ftX, ftY; case = :middle);
Equipercentile(ftX, ftY; case = :middle).table.eYx |> print
ansEp = [0.3742266165402808, 2.8668002442948866, 4.377679952675711, 5.589899901936712, 6.5719466795280965, 7.675963305503476, 8.779402720766932, 9.821396984022385, 10.832444129170566, 11.807020939106527, 12.762424501360979, 13.674737884179626, 14.606852126241565, 15.531810687908973, 16.384026327429506, 17.269396166667576, 18.13458065086816, 18.993135811329992, 19.859585529606615, 20.69055191799923, 21.57870500928172, 22.426825738450507, 23.24982055643067, 24.084693805207607, 24.971262459510967, 25.86435538130699, 26.775635117977668, 27.793162605037413, 28.849521797543677, 29.862255991592047, 30.864456663583688, 31.869824776800975, 32.91214157149059, 33.951819244636646, 34.89830208873431, 35.799993897278775, 36.82494259856727, 37.87655236332893, 38.90275680579959, 40.018418515159]
all(resEp.table.eYx .≈ ansEp)
# Linear Equating
reslin = Linear(ftX, ftY)
coef(reslin)
# Mean Equating
Mean(ftX, ftY) |> coef

# Run equating under the SG design
using CSVFiles, DataFrames, RCall
@rimport equate as requate
ACTmath = DataFrame!(load("data/ACTmath.csv"))
X = fill.(ACTmath.scale, ACTmath.xcount) |> Iterators.flatten |> collect
Y = fill.(ACTmath.scale, ACTmath.ycount) |> Iterators.flatten |> collect
ftX = freqtab(X)
ftY = freqtab(Y)
rftX = requate.freqtab(X)
rftY = requate.freqtab(Y)

# equipercentile
resE = Equipercentile(ftX, ftY; case = :lower)
resE.table
resEr = requate.equate(rftX, rftY, type = "e")
resEr["concordance"]
resEr["concordance"][2] .- resE.table.eYx
# linear
resL = Linear(ftX, ftY)
resLr = requate.equate(rftX, rftY, type = "l")
resLr["concordance"][2] .- resL.table[:, r"Yx"][:, 1]
# Mean
resM = Mean(ftX, ftY)
resMr = requate.equate(rftX, rftY, type = "m")
resMr["concordance"][2] .- resM.table[:, r"Yx"][:, 1]\

# Log Linear Smoothing
using CSVFiles, DataFrames, GLM
KBneatX = DataFrame!(load("data/KBneatX.csv"))
ftX = freqtab(KBneatX.total)
ftXsmoothed, fit1 = presmoothing(ftX, LogLinearFormula(6))
predict(fit1)
using Plots
plot(ftXsmoothed.tab.scale, predict(fit1), label = "Smoothed degree = 6")
plot!(ftX.tab.scale, ftX.tab.freq; label = "observed probability")
# link[https://github.com/JuliaStats/StatsModels.jl/blob/master/docs/src/formula.md]

# Kernel Smoothing
using CSVFiles, DataFrames
ACTmath = DataFrame!(load("data/ACTmath.csv"))
X = ExpandTable(ACTmath.scale, ACTmath.xcount)
ftX = freqtab(X; scale = 0:1:40)
smftX = presmoothing(ftX, LogLinearFormula(4))
# Choose bandwidth
optimalbwidth = EstBandwidth(smftX)
optimalbwidth.minimizer[1]

KftX = KernelSmoothing(smftX; hX = exp(optimalbwidth.minimizer[1]))
KftX = KernelSmoothing(smftX; hX = 100)
KftX = KernelSmoothing(smftX; hX = 0.33)

using Plots
plot(KftX)

# Smoothed SG
using CSVFiles, DataFrames
ACTmath = DataFrame!(load("data/ACTmath.csv"))
X = ExpandTable(ACTmath.scale, ACTmath.xcount)
Y = ExpandTable(ACTmath.scale, ACTmath.ycount)
ftX = freqtab(X)
ftY = freqtab(Y)
smX = presmoothing(ftX, LogLinearFormula(4));
@time fml = @LogLinearFormula(4);
@time smY = presmoothing(ftY, @LogLinearFormula(4));
@time smY = presmoothing(ftY, fml);
Equipercentile(ftX, ftY)

presmoothing(ftY, 6)

using GLM
deviance(smY.fit)
dof_residual(smY.fit)
using StatsBase
aic(smY.fit)
bic(smY.fit)
aicc(smY.fit)
loglikelihood(smY.fit)
deviance(smY.fit)

ftX = KernelSmoothing(freqtab(X))
ftY = KernelSmoothing(freqtab(Y))
Equipercentile(ftX, ftY)

# Nonequivalent group design
using CSVFiles, DataFrames
KBneatX = DataFrame!(load("test/data/KBneatX.csv"))
KBneatY = DataFrame!(load("test/data/KBneatY.csv"))

ftX = freqtab(KBneatX.total, KBneatX.anchor)
ftY = freqtab(KBneatY.total, KBneatY.anchor)
sum(ftX.marginal)
heatmap(ftX.marginal, color = :plasma)

# Run equiting functions under the NEAT design and results comparison.
using RCall
@rimport equate as requate
rftX = requate.freqtab(DataFrame(total = KBneatX.total, anchor = KBneatX.anchor))
rftY = requate.freqtab(DataFrame(total = KBneatY.total, anchor = KBneatY.anchor))
#Tucker # Matched
resTk = Tucker(ftX, ftY)
resTkr = requate.equate(rftX, rftY, type = "l", method = "tucker")
plot(resTk.table.scaleX, resTk.table.lYx; label = "Tucker", xlabel = "scale X", ylabel = "scale Y")

# Levine # common = :internal is matched
resLv = LevineCongeneric(ftX, ftY; common = :internal)
resTkr = requate.equate(rftX, rftY, type = "l", method = "levine")
# external
resLv = LevineCongeneric(ftX, ftY; common = :external)
resTkr = requate.equate(rftX, rftY, type = "l", method = "levine", internal = false)

# Chained Linear # Matched
ObservableStats(ftY)
resCL = ChainedLinear(ftX, ftY)
resCLr = requate.equate(rftX, rftY, type = "l", method = "chained")
plot!(resCL.table.scaleX, resCL.table.lYx; label = "Chained Linear", xlabel = "scale X", ylabel = "scale Y")

# Braun & Holland # Not matched
resBH = BraunHolland(ftX, ftY)
resBHr = requate.equate(rftX, rftY, type = "l", method = "braun/holland")
plot!(resBH.table.scaleX, resBH.table.lYx; label = "Braun Holland", xlabel = "scale X", ylabel = "scale Y")

# Frequency estimation
ftm = ftX.marginal
resFE = FrequencyEstimation(ftX, ftY; case = :lower)
resFEr = requate.equate(rftX, rftY, type = "e", method = "frequency estimation")
resFEr["concordance"][2] .- resFE.table.eYx
plot!(resFE.table.eYx, resFE.table.scaleY; label = "Frequency Estimation", xlabel = "scale X", ylabel = "scale Y")

# Chained Equipercentile
resCE = ChainedEquipercentile(ftX, ftY; case = :upper)
resCEr = requate.equate(rftX, rftY, type = "e", method = "chained")
resCEr["concordance"][2] .- resCE.table.eYx
plot!(resCE.table.eYx, resCE.table.scaleX; label = "Chained Equipercentile", xlabel = "scale X", ylabel = "scale Y", legend = :topleft)

# smoothed NEAT
# @profview 
ft1 = presmoothing(ftX, LogLinearFormula(6), LogLinearFormula(6))
ft2 = presmoothing(ftY, LogLinearFormula(6), LogLinearFormula(6))
ChainedLinear(ft1, ft2)

# Extract coefficient
using CSVFiles, DataFrames
KBneatX = DataFrame!(load("data/KBneatX.csv"))
KBneatY = DataFrame!(load("data/KBneatY.csv"))

ftX = freqtab(KBneatX.total, KBneatX.anchor)
ftY = freqtab(KBneatY.total, KBneatY.anchor)

resTk = Tucker(ftX, ftY)
coef(resTk)

resFE = FrequencyEstimation(ftX, ftY)

resCE = ChainedEquipercentile(ftX, ftY)
coef(resCE)

# Plot recipes
using Plots
using CSVFiles, DataFrames
KBneatX = DataFrame!(load("data/KBneatX.csv"))
KBneatY = DataFrame!(load("data/KBneatY.csv"))

ftX = freqtab(KBneatX.total, KBneatX.anchor)
ftY = freqtab(KBneatY.total, KBneatY.anchor)

resTk = Tucker(ftX, ftY)
plot(resTk)

ftX = freqtab(KBneatX.total)
KftX = KernelSmoothing(ftX; hX = 0.66)
plot(KftX)

smftX = presmoothing(ftX, LogLinearFormula(6))
plot(smftX)
plot(smftX, smftX.fit)

plot(ftX)

# Summary Stats
using CSVFiles, DataFrames
KBneatX = DataFrame!(load("data/KBneatX.csv"))
KBneatY = DataFrame!(load("data/KBneatY.csv"))

ftX = freqtab(KBneatX.total, KBneatX.anchor)
ftY = freqtab(KBneatY.total, KBneatY.anchor)

test = SummaryStats(ftX);

# Presmooth model comparison
ftX = freqtab(KBneatX.total)

using CSV, Plots
ACTmath = DataFrame!(CSV.File("test/data/ACTmath.csv"))
X = ExpandTable(ACTmath.scale, ACTmath.xcount)
Y = expandtable(ACTmath.scale, ACTmath.ycount)
ftX = freqtab(X; scale = 0:1:40)
ftY = freqtab(Y; scale = 0:1:40)

plot()

Equipercentile(ftX, ftY)
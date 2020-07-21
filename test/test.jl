# Equipercentile Equating
using CSVFiles, DataFrames
ACTmath = DataFrame!(load("data/ACTmath.csv"))
X = fill.(ACTmath.scale, ACTmath.xcount) |> Iterators.flatten |> collect
Y = fill.(ACTmath.scale, ACTmath.ycount) |> Iterators.flatten |> collect
ftX = freqtab(X); ftY = freqtab(Y);
Equipercentile(ftX, ftY; case = :upper) |> print
Equipercentile(ftX, ftY; case = :lower) |> print
Equipercentile(ftX, ftY; case = :both) |> print
resEp = Equipercentile(ftX, ftY; case = :middle);
Equipercentile(ftX, ftY; case = :middle).table.eYx |> print
ansEp = [0.3742266165402808, 2.8668002442948866, 4.377679952675711, 5.589899901936712, 6.5719466795280965, 7.675963305503476, 8.779402720766932, 9.821396984022385, 10.832444129170566, 11.807020939106527, 12.762424501360979, 13.674737884179626, 14.606852126241565, 15.531810687908973, 16.384026327429506, 17.269396166667576, 18.13458065086816, 18.993135811329992, 19.859585529606615, 20.69055191799923, 21.57870500928172, 22.426825738450507, 23.24982055643067, 24.084693805207607, 24.971262459510967, 25.86435538130699, 26.775635117977668, 27.793162605037413, 28.849521797543677, 29.862255991592047, 30.864456663583688, 31.869824776800975, 32.91214157149059, 33.951819244636646, 34.89830208873431, 35.799993897278775, 36.82494259856727, 37.87655236332893, 38.90275680579959, 40.018418515159]
all(resEp.table.eYx .≈ ansEp)
# Linear Equating
Linear(ftX, ftY).table.lYx |> print

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
smftX = presmoothing(ftX, LogLinearFormula(6))
# Choose bandwidth
optimalbwidth = EstBandwidth(smftX)
optimalbwidth.minimizer[1]

KftX = KernelSmoothing(smftX; hX = exp(optimalbwidth.minimizer[1]))
KftX = KernelSmoothing(smftX; hX = 100)
KftX = KernelSmoothing(smftX; hX = 0.33)

using Plots
plot(KftX)

# Smoothed SG
using CSV
ACTmath = CSV.read("test/ACTmath.csv")
X = fill.(ACTmath.scale, ACTmath.xcount) |> Iterators.flatten |> collect
Y = fill.(ACTmath.scale, ACTmath.ycount) |> Iterators.flatten |> collect
ftX, fit1 = presmoothing(freqtab(X); fml = LogLinearFormula(4))
ftY, fit2 = presmoothing(freqtab(Y); fml = LogLinearFormula(4))
Equipercentile(ftX, ftY)

ftX = KernelSmoothing(freqtab(X))
ftY = KernelSmoothing(freqtab(Y))
Equipercentile(ftX, ftY)

# Nonequivalent group design
using CSVFiles, DataFrames
KBneatX = DataFrame!(load(("test/KBneatX.csv"))
KBneatY = DataFrame!(load(("test/KBneatY.csv"))

ftX = freqtab(KBneatX.total, KBneatX.anchor)
ftY = freqtab(KBneatY.total, KBneatY.anchor)
sum(ftX.marginal)
heatmap(ftX.marginal, color = :plasma)

#Tucker
resTk = Tucker(ftX, ftY)
plot(resTk.table.scaleX, resTk.table.lYx; label = "Tucker", xlabel = "scale X", ylabel = "scale Y")

# Chained Linear
ObservableStats(ftY)
resCL = ChainedLinear(ftX, ftY)
plot!(resCL.table.scaleX, resCL.table.lYx; label = "Chained Linear", xlabel = "scale X", ylabel = "scale Y")

# Frequency estimation
ftm = ftX.marginal
resFE = FrequencyEstimation(ftX, ftY)
plot!(resFE.table.eYx, resFE.table.scaleY; label = "Frequency Estimation", xlabel = "scale X", ylabel = "scale Y")

# Braun & Holland
resBH = BraunHolland(ftX, ftY)
plot!(resBH.table.scaleX, resBH.table.lYx; label = "Braun Holland", xlabel = "scale X", ylabel = "scale Y")

# Chained Equipercentile
resCE = ChainedEquipercentile(freqtab(KBneatX.total, KBneatX.anchor), freqtab(KBneatY.total, KBneatY.anchor))
plot!(resCE.table.eYx, resCE.table.scaleY; label = "Chained Equipercentile", xlabel = "scale X", ylabel = "scale Y", legend = :topleft)

# smoothed NEAT
ftX, fit1 = presmoothing(freqtab(KBneatX.total); fml = LogLinearFormula(6))
ftXV, fit2 = presmoothing(freqtab(KBneatX.anchor); fml = LogLinearFormula(6))
ftY, fit3 = presmoothing(freqtab(KBneatY.total); fml = LogLinearFormula(6))
ftYV, fit4 = presmoothing(freqtab(KBneatY.anchor); fml = LogLinearFormula(6))

ft1 = freqtab(ftX, ftXV)
ft2 = freqtab(ftY, ftYV)
ChainedEquipercentile(ft1, ft2)

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
plot(smftX, smfit)
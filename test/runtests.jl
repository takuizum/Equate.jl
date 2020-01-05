using Equate
using Test
using CSV
dir = replace(pwd(), "test" => "")
ACTmath = CSV.read(dir*"data/ACTmath.csv")
KBneatX = CSV.read(dir*"data/KBneatX.csv")
KBneatY = CSV.read(dir*"data/KBneatY.csv")
# using RCall
# R"""
# install.packages("equate")
# library(equate)
# KBneatX = KBneat$x
# KBneatY = KBneat$y
# """
# @rget ACTmath
# @rget KBneatX
# @rget KBneatY

# answer vector
ansEp = [0.3742266165402808, 2.8668002442948866, 4.377679952675711, 5.589899901936712, 6.5719466795280965, 7.675963305503476, 8.779402720766932, 9.821396984022385, 10.832444129170566, 11.807020939106527, 12.762424501360979, 13.674737884179626, 14.606852126241565, 15.531810687908973, 16.384026327429506, 17.269396166667576, 18.13458065086816, 18.993135811329992, 19.859585529606615, 20.69055191799923, 21.57870500928172, 22.426825738450507, 23.24982055643067, 24.084693805207607, 24.971262459510967, 25.86435538130699, 26.775635117977668, 27.793162605037413, 28.849521797543677, 29.862255991592047, 30.864456663583688, 31.869824776800975, 32.91214157149059, 33.951819244636646, 34.89830208873431, 35.799993897278775, 36.82494259856727, 37.87655236332893, 38.90275680579959, 40.018418515159]
ansLin = [-1.5433492744496915, -0.45472777167758416, 0.6338937310945232, 1.7225152338666305, 2.811136736638738, 3.899758239410845, 4.9883797421829525, 6.07700124495506, 7.165622747727166, 8.254244250499275, 9.342865753271383, 10.43148725604349, 11.520108758815596, 12.608730261587704, 13.697351764359812, 14.785973267131919, 15.874594769904025, 16.96321627267613, 18.05183777544824, 19.140459278220348, 20.229080780992454, 21.317702283764564, 22.40632378653667, 23.494945289308777, 24.583566792080884, 25.67218829485299, 26.7608097976251, 27.849431300397207, 28.938052803169313, 30.026674305941423, 31.11529580871353, 32.20391731148564, 33.29253881425774, 34.38116031702985, 35.46978181980195, 36.558403322574065, 37.64702482534618, 38.73564632811828, 39.82426783089039, 40.91288933366249]

@testset "Equate.jl" begin
    # SG design
    X = fill.(ACTmath.scale, ACTmath.xcount) |> Iterators.flatten |> collect
    Y = fill.(ACTmath.scale, ACTmath.ycount) |> Iterators.flatten |> collect
    ftX = freqtab(X); ftY = freqtab(Y)
    resEp = Equipercentile(ftX, ftY; case = :middle)
    @test all(resEp.table.eYx .≈ ansEp)

    resLin = Linear(ftX, ftY)
    @test all(resLin.table.lYx .≈ ansLin)
    # Smoothing
    ftXsmoothed, fit1 = presmoothing(ftX; fml = LogLinearFormula(6))
    BandwidthOpt = EstBandwidth(ftX)
    ftXKsmoothed = KernelSmoothing(ftX; hX = 0.622)

    # NEAT design |>
    ftX = freqtab(KBneatX.total, KBneatX.anchor)
    ftY = freqtab(KBneatY.total, KBneatY.anchor)
    Tucker(ftX, ftY)
    BraunHolland(ftX, ftY)
    ChainedLinear(ftX, ftY)
    ChainedEquipercentile(ftX, ftY)

    # Smoothed SG
    # Presmoothed Equipercentile
    ftX, fit1 = presmoothing(freqtab(X); fml = LogLinearFormula(4))
    ftY, fit2 = presmoothing(freqtab(Y); fml = LogLinearFormula(4))
    Equipercentile(ftX, ftY)
    # Kernel Equating
    ftX = KernelSmoothing(freqtab(X))
    ftY = KernelSmoothing(freqtab(Y))
    Equipercentile(ftX, ftY)

    # Smoothed NEAT
    ftX = freqtab(KBneatX.total, KBneatX.anchor)
    ftY = freqtab(KBneatY.total, KBneatY.anchor)
    resTk = Tucker(ftX, ftY)
    resCL = ChainedLinear(ftX, ftY)
    resFE = FrequencyEstimation(ftX, ftY)
    resCE = ChainedEquipercentile(ftX, ftY)
    resBH = BraunHolland(ftX, ftY)

end

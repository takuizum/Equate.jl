using Equate
using Test

@testset "Equate.jl" begin
    using RCall
    R"""
    # install.packages("equate")
    library(equate)
    """
    @rget ACTmath
    @rget KBneat
    #
    # using CSV
    # ACTmath = CSV.read("ACTmath.csv")
    # KBneatX = CSV.read("KBneatX.csv")
    # KBneatY = CSV.read("KBneatY.csv")

    # SG design
    X = fill.(ACTmath.scale, ACTmath.xcount) |> Iterators.flatten |> collect
    Y = fill.(ACTmath.scale, ACTmath.ycount) |> Iterators.flatten |> collect
    ftX = freqtab(X); ftY = freqtab(Y)
    Equipercentile(ftX, ftY; case = :middle)
    Linear(ftX, ftY)
    # Smoothing
    ftXsmoothed, fit1 = presmoothing(ftX; fml = LogLinearFormula(6))
    BandwidthOpt = EstBandwidth(ftX)
    ftXKsmoothed = KernelSmoothing(ftX; hX = 0.622)

    # NEAT design
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

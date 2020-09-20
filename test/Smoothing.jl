using StatsBase

@testset "presmoothing" begin
    # df = 1
    fml1 = LogLinearFormula(1)
    smtab = presmoothing(sgtabX, fml1)
    m = smtab.tab.scale'smtab.tab.prob
    @test m ≈ mean(smtab.raw)
    # df = 2
    fml2 = LogLinearFormula(2)
    smtab = presmoothing(sgtabX, fml2)
    m = smtab.tab.scale'smtab.tab.prob
    s = ((smtab.tab.scale .- m).^2)'smtab.tab.prob
    @test s ≈ var(smtab.raw; corrected = false)
    # df = 3
    fml3 = LogLinearFormula(3)
    smtab = presmoothing(sgtabX, fml3)
    m = smtab.tab.scale'smtab.tab.prob
    s = ((smtab.tab.scale .- m).^2)'smtab.tab.prob
    sk = ((smtab.tab.scale .- m).^3)'smtab.tab.prob / sqrt(s)^3
    # @test sk ≈ skewness(smtab.raw)
    @test round(sk; digits = 5) ≈ round(skewness(smtab.raw); digits = 5)
    # df = 4
    fml4 = LogLinearFormula(4)
    smtab = presmoothing(sgtabX, fml4)
    m = smtab.tab.scale'smtab.tab.prob
    s = ((smtab.tab.scale .- m).^2)'smtab.tab.prob
    ku = ((smtab.tab.scale .- m).^4)'smtab.tab.prob / s^2 -3
    @test ku ≈ kurtosis(smtab.raw)
end

@testset "KernelSmoothing (Gaussian)" begin
    fml = LogLinearFormula(6)
    smtab = presmoothing(sgtabX, fml)
    bbw = EstBandwidth(smtab)
    @test bbw.ls_success
    ktab = KernelSmoothing(smtab; hX = exp(bbw.minimizer[1]))
    @test typeof(ktab) == Equate.KernelFreqTab
end


using StatsBase

@testset "presmoothing" begin
    # df = 1
    smtab = presmoothing(sgtabX, Equate.fml₁)
    m = smtab.tab.scale'smtab.tab.prob
    @test m ≈ mean(smtab.raw)
    # df = 2
    smtab = presmoothing(sgtabX, Equate.fml₂)
    m = smtab.tab.scale'smtab.tab.prob
    s = ((smtab.tab.scale .- m).^2)'smtab.tab.prob
    @test s ≈ var(smtab.raw; corrected = false)
    # df = 3
    smtab = presmoothing(sgtabX, Equate.fml₃)
    m = smtab.tab.scale'smtab.tab.prob
    s = ((smtab.tab.scale .- m).^2)'smtab.tab.prob
    sk = ((smtab.tab.scale .- m).^3)'smtab.tab.prob / sqrt(s)^3
    # @test sk ≈ skewness(smtab.raw)
    @test round(sk; digits = 5) ≈ round(skewness(smtab.raw); digits = 5)
    # df = 4
    smtab = presmoothing(sgtabX, Equate.fml₄)
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

@testset "Examine smoothing" begin
    fit = presmoothing(sgtabX, 6)
    @test fit isa DataFrame
end


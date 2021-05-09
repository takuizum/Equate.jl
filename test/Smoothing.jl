using StatsBase

@testset "presmoothing" begin
    # df = 1
    smtab = presmoothing(sgtabX, Equate.fml₁)
    m = smtab.table.scale'smtab.table.prob
    @test m ≈ mean(smtab.raw)
    # df = 2
    smtab = presmoothing(sgtabX, Equate.fml₂)
    m = smtab.table.scale'smtab.table.prob
    s = ((smtab.table.scale .- m).^2)'smtab.table.prob
    @test s ≈ var(smtab.raw; corrected = false)
    # df = 3
    smtab = presmoothing(sgtabX, Equate.fml₃)
    m = smtab.table.scale'smtab.table.prob
    s = ((smtab.table.scale .- m).^2)'smtab.table.prob
    sk = ((smtab.table.scale .- m).^3)'smtab.table.prob / sqrt(s)^3
    # @test sk ≈ skewness(smtab.raw)
    @test round(sk; digits = 5) ≈ round(skewness(smtab.raw); digits = 5)
    # df = 4
    smtab = presmoothing(sgtabX, Equate.fml₄)
    m = smtab.table.scale'smtab.table.prob
    s = ((smtab.table.scale .- m).^2)'smtab.table.prob
    ku = ((smtab.table.scale .- m).^4)'smtab.table.prob / s^2 -3
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


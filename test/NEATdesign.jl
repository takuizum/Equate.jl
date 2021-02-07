using Equate
using Statistics, StatsBase

# methods = [
#     "Tucker"
#     "LevineCongeneric"
#     "ChainedLinear"
#     "FrequencyEstimation"
#     "BraunHolland"
#     "ChainedEquipercentile"
# ]

@testset "Tucker" begin
    res = Tucker(neattabX, neattabY)
    @test round(coef(res)[1]; digits = 4) == 1.0292
    @test round(coef(res)[2]; digits = 4) == .5378
    # popXonY = ExpandTable(res.table.lYx, neattabX.tabX.freq)
    # @test round(mean(popXonY); digits = 4) == 16.8153
end

@testset "LevineCongeneric (external) " begin
    res = LevineCongeneric(neattabX, neattabY; common = :external)
    @test round(coef(res)[1]; digits = 4) == 1.0168
    @test round(coef(res)[2]; digits = 4) == .3257
    # popXonY = ExpandTable(res.table.lYx, neattabX.tabX.freq)
    # @test round(mean(popXonY); digits = 4) == 16.2485
end

@testset "LevineCongeneric (internal) " begin
    res = LevineCongeneric(neattabX, neattabY; common = :internal)
    @test round(coef(res)[1]; digits = 4) == 1.0110
    @test round(coef(res)[2]; digits = 4) == .2514
end

@testset "ChainedLinear" begin
    res = ChainedLinear(neattabX, neattabY)
    @test round(coef(res)[1]; digits = 4) == 1.0213
    @test round(coef(res)[2]; digits = 4) == .3937
end

# Invoke the percentile rank methods
# The results in Kolen & Brennan seems to be wrong.
@testset "BraunHolland (R-equate)" begin
    res = BraunHolland(neattabX, neattabY)
    @test round(coef(res)[1]; digits = 4) == 1.0067
    @test round(coef(res)[2]; digits = 4) == 0.8976
    # popXonY = ExpandTable(res.table.lYx, neattabX.tabX.freq)
    # @test round(mean(popXonY); digits = 4) == 16.8329
    # @test round(std(popXonY; corrected = false); digits = 4) == 6.6017
    # @test round(skewness(popXonY); digits = 4) == .5799
end

# @testset "BraunHolland (Kolen & Brennan)" begin
#     res = BraunHolland(dummytabX, dummytabY)
#     @test round(coef(res)[1]; digits = 4) == 1.0331
#     @test round(coef(res)[2]; digits = 4) == -.1927
#     popXonY = ExpandTable(res.table.lYx, Int64.(dummytabX.tabX.freq))
#     @test round(mean(popXonY); digits = 4) == 16.8329
#     @test round(std(popXonY; corrected = false); digits = 4) == 6.6017
#     @test round(skewness(popXonY); digits = 4) == .5799
# end

@testset "FrequencyEstimation" begin
    res = FrequencyEstimation(dummytabX, dummytabY; case = :upper, w₁ = 1.0)
    @test round.(res.table.eYx; digits = 2) == [-0.02, 0.83, 1.76, 2.92, 3.98, 5.00]
end

@testset "ChainedEquipercentile" begin
    res = ChainedEquipercentile(neattabX, neattabY; case = :upper)
    # Largest 3 values in the concordance table are matched.
    @test round.(res.table.eYx[end-2:end]; digits = 4) == [34.3125, 35.4125, 36.0938]
end


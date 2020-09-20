using Equate
using Statistics

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

@testset "LevineCongeneric (internal) " begin
    res = LevineCongeneric(neattabX, neattabY; common = :internal)
    @test round(coef(res)[1]; digits = 4) == 1.0110
    @test round(coef(res)[2]; digits = 4) == .2514
    # popXonY = ExpandTable(res.table.lYx, neattabX.tabX.freq)
    # @test round(mean(popXonY); digits = 4) == 16.2485
end

@testset "ChainedLinear" begin
    res = ChainedLinear(neattabX, neattabY)
    @test round(coef(res)[1]; digits = 4) == 1.0213
    @test round(coef(res)[2]; digits = 4) == .3937
end

# Invoke the percentile rank methods
# The results in Kolen & Brennan seems to be wrong.
# @testset "BraunHolland" begin
#     res = BraunHolland(neattabX, neattabY)
#     # @test round(coef(res)[1]; digits = 4) == 1.0213
#     # @test round(coef(res)[2]; digits = 4) == .3937
#     popXonY = ExpandTable(res.table.lYx, neattabX.tabX.freq)
#     @test round(mean(popXonY); digits = 4) == 16.8329
#     @test round(std(popXonY; corrected = false); digits = 4) == 6.6017
#     @test round(skewness(popXonY); digits = 4) == .5799
# end

# @testset "BraunHolland" begin
#     res = BraunHolland(dummytabX, dummytabY)
#     @test round(coef(res)[1]; digits = 4) == 1.0331
#     @test round(coef(res)[2]; digits = 4) == -.1927
# end

# @testset "FrequencyEstimation" begin
#     res = FrequencyEstimation(dummytabX, dummytabY)
# end


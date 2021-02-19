include("testdata.jl")

@testset "SGcoef" begin
    c = coef(Linear(sgtabX, sgtabY))
    @test c.slope ≈  1.08862150
    @test c.intercept ≈ -2.63197077
end

@testset "NEATcoef" begin
    c = coef(ChainedLinear(neattabX, neattabY))
    @test c.slope ≈ 1.0212716885774726
    @test c.intercept ≈ 0.39367983874242185
end
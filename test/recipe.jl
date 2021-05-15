# plot.jl

recipes_testing = false

if recipes_testing
    using StatsPlots, Equate

    include("readdata.jl")
    sgtabX

    plot(sgtabX)
end

@testset "recipes" begin
    println("Recipes testing is a computer intensive. So skipped.")
    @test true
end
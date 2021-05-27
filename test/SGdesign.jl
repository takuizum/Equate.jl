using Equate
using Statistics

@testset "Mean" begin
    res = Mean(sgtabX, sgtabY)
    popXonY = ExpandTable(res.table.lYx, ACTmath.xcount)
    @test round(mean(popXonY); digits = 4) == 18.9798
end

@testset "Linear" begin
    res = Linear(sgtabX, sgtabY)
    popXonY = ExpandTable(res.table.lYx, ACTmath.xcount)
    @test round(mean(popXonY); digits = 4) == 18.9798
    @test round(std(popXonY; corrected = false); digits = 4) == 8.9394 # 8.9393 in KB p.50
end

@testset "Equipercentile (lower)" begin
    res = Equipercentile(sgtabX, sgtabY)
    popXonY = ExpandTable(res.table.eYx, ACTmath.xcount)
    @test round(mean(popXonY); digits = 4) == 18.9799
    @test round(std(popXonY; corrected = false); digits = 4) == 8.9352
end

@testset "Percentile Rank Function" begin
    P̂ = round2.(Equate.PRF.(0:1:40, Ref(sgtabX)); digits = 2)
    @test sum(KBp48 .== P̂) == 41
end


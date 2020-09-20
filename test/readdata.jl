using Test
using Equate
using CSV, Statistics, DataFrames

ACTmath = DataFrame!(CSV.File("data/ACTmath.csv"))
KBneatX = DataFrame!(CSV.File("data/KBneatX.csv"))
KBneatY = DataFrame!(CSV.File("data/KBneatY.csv"))

ACTmathX = ExpandTable(ACTmath.scale, ACTmath.xcount)
ACTmathY = ExpandTable(ACTmath.scale, ACTmath.ycount)

@testset "ExpandTable" begin
    @test ExpandTable([0, 1, 2], [2, 2, 2]) == [0, 0, 1, 1, 2, 2]
    @test expandtable([0, 1, 2], [2, 2, 2]) == [0, 0, 1, 1, 2, 2]
end

@testset "DataMoments" begin
    @test round(mean(ACTmathX); digits = 4) == 19.8524
    @test round(std(ACTmathX; corrected = false); digits = 4) == 8.2116
    @test round(mean(ACTmathY); digits = 4) == 18.9798
    @test round(std(ACTmathY; corrected = false); digits = 4) == 8.9393
end

sgtabX = freqtab(ACTmathX; scale = 0:1:40)
sgtabY = freqtab(ACTmathY; scale = 0:1:40)

neattabX = freqtab(KBneatX.total, KBneatX.anchor; scaleX = 0:1:36, scaleV = 0:1:12)
neattabY = freqtab(KBneatY.total, KBneatY.anchor; scaleX = 0:1:36, scaleV = 0:1:12)

# Dummy table for NEAT equipercentile method
dummytabX = Equate.SGFreqTab(
    DataFrame(scale = [0, 1, 2, 3, 4, 5], freq = [.1, .15, .25, .25, .15, .10] .* 100, prob = [.1, .15, .25, .25, .15, .10], cumprob = cumsum([.1, .15, .25, .25, .15, .10]) ), 
    DataFrame(scale = [0, 1, 2, 3], freq = [.2, .4, .2, .2] .* 100, prob = [.2, .4, .2, .2], cumprob = cumsum([.2, .4, .2, .2])), 
    ExpandTable([0, 1, 2, 3, 4, 5], Int64.([.1, .15, .25, .25, .15, .10] .* 100)), 
    ExpandTable([0, 1, 2, 3], Int64.([.2, .4, .2, .2] .* 100)), 
    1, 
    1, 
    [.04 .04 .02 .00
     .04 .08 .02 .01
     .06 .12 .05 .02
     .03 .12 .05 .05
     .02 .03 .04 .06
     .01 .01 .02 .06] .* 100
)

dummytabY = Equate.SGFreqTab(
    DataFrame(scale = [0, 1, 2, 3, 4, 5], freq = [.08, .2, .22, .25, .15, .1] .* 100, prob = [.08, .2, .22, .25, .15, .1], cumprob = cumsum([.08, .2, .22, .25, .15, .1]) ), 
    DataFrame(scale = [0, 1, 2, 3], freq = [.2, .2, .4, .2] .* 100, prob = [.2, .2, .4, .2], cumprob = cumsum([.2, .2, .4, .2])), 
    ExpandTable([0, 1, 2, 3, 4, 5], Int64.([.08, .2, .22, .25, .15, .1] .* 100)), 
    ExpandTable([0, 1, 2, 3], Int64.([.2, .2, .4, .2] .* 100)), 
    1, 
    1, 
    [.04 .03 .01 .00
     .07 .05 .07 .01
     .03 .05 .12 .02
     .03 .04 .13 .05
     .02 .02 .05 .06
     .01 .01 .02 .06] .* 100
)


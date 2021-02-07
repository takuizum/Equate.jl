using Test
using Equate
using CSV, Statistics, DataFrames

# ACTmath = DataFrame!(CSV.File("test/data/ACTmath.csv"))
ACTmath = DataFrame([0 0 0
1 1 1
2 1 3
3 3 13
4 9 42
5 18 59
6 59 95
7 67 131
8 91 158
9 144 161
10 149 194
11 192 164
12 192 166
13 192 197
14 201 177
15 204 158
16 217 169
17 181 132
18 184 158
19 170 151
20 201 134
21 147 137
22 163 122
23 147 110
24 140 116
25 147 132
26 126 104
27 113 104
28 100 114
29 106 97
30 107 107
31 91 88
32 83 80
33 73 79
34 72 70
35 75 61
36 50 48
37 37 47
38 38 29
39 23 32
40 15 12], ["scale","xcount","ycount"])

KBneatX = DataFrame!(CSV.File("test/data/KBneatX.csv"))
KBneatY = DataFrame!(CSV.File("test/data/KBneatY.csv"))

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


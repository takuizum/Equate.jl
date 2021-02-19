include("testdata.jl")

@testset "expandtable" begin
    dat = expandtable(ACTmath.scale, ACTmath.xcount)
    tab = freqtab(dat)
    @test tab.tab.freq == ACTmath.xcount[begin+1:end]
end

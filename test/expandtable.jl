include("testdata.jl")

@testset "expandtable(Single)" begin
    dat = expandtable(ACTmath.scale, ACTmath.xcount)
    ft = freqtab(dat)
    @test ft.tab.freq == ACTmath.xcount[begin+1:end]
end

include("testdata.jl")

@testset "expandtable(Single)" begin
    dat = expandtable(ACTmath.scale, ACTmath.xcount)
    ft = freqtab(dat)
    @test ft.table.freq == ACTmath.xcount[begin+1:end]
end

using Bootstrap

@testset "recalcuate!" begin
    # SG
    res = Equipercentile(sgtabX, sgtabY)
    res.data.X.raw = collect(1:40)
    recalculate!(res.data.X)
    @test res.data.X.stats.μ == 20.5
    @test res.data.X.table.freq == map(j -> count(i -> i == j, collect(1:40)), res.data.X.table.scale)

    # NEAT
    test = copy(neattabX)
    test.rawX = collect(0:10)
    recalculate!(test)
    @test test.tableX != neattabX.tableX
    @test test.tableV == neattabX.tableV
end

@testset "Base.copy for FreqTab" begin
    res = Linear(sgtabX, sgtabY)
    test = copy(res.data)
    merge(test, (X = freqtab([1,2,3,4,5]), ))
    @test res.data.X != test.X
    test.X.raw = nothing
    @test res.data.X.raw != test.X.raw
end

@testset "Base.copy for SGEquateResult" begin
    res = Linear(sgtabX, sgtabY)
    test = copy(res)
    initialize!(test.data.X)
    @test test.data.X != res.data.X
end

@testset "Base.copy for NEATFreqTab" begin
    test = copy(neattabX)    
    test.statsX = 0
    @test neattabX.statsX != test.statsX
end

@testset "Base.copy for NEATEquateResult" begin
    res = BraunHolland(neattabX, neattabY)
    test = copy(res)
    @test true # Not good!
end

@testset "Bootstrap.draw!" begin
    seed!(1234)
    res = Linear(sgtabX, sgtabY)
    o = copy(res.data)
    Bootstrap.draw!(res.data, o)
    @test res.data.X.stats != o.X.stats
    @test res.data.Y.stats != o.Y.stats
end

@testset "Bootstrap.bootstrap for SG" begin
    seed!(1234)
    res = Linear(sgtabX, sgtabY)
    bsres = bootstrap(x -> coef(Linear(x...)), res.data, BasicSampling(10))
    @test size(Bootstrap.estimate_summary(bsres)) == (2, 3)
    seed!(1234)
    res = Mean(sgtabX, sgtabY)
    bsres = bootstrap(x -> coef(Mean(x...)), res.data, BasicSampling(10))
    @test size(Bootstrap.estimate_summary(bsres)) == (2, 3)
    seed!(1234)
    res = Equipercentile(sgtabX, sgtabY)
    bsres = bootstrap(x -> coef(Equipercentile(x...)), res.data, BasicSampling(10))
    @test size(Bootstrap.estimate_summary(bsres)) == (41, 3)
end

@testset "Bootstrap.bootstrap for NEAT" begin
    for f in [Tucker, LevineObservedScore, LevineTrueScore, ChainedLinear, BraunHolland, FrequencyEstimation, ChainedEquipercentile]
        res = f(neattabX, neattabY)
        bsres = bootstrap(x -> coef(f(x...)), res.data, BasicSampling(50))
        @test true
    end
end


# Other bootstrap method

# res = Linear(sgtabX, sgtabY)
# bootstrap(x -> coef(Linear(x...)), res.data, BasicSampling(10)) # OK
# bootstrap(x -> coef(Linear(x...)), res.data, AntitheticSampling(10)) # Only for vectored data.
# bootstrap(x -> coef(Linear(x...)), res.data, BalancedSampling(10)) #
# bootstrap(x -> coef(Linear(x...)), res.data, MaximumEntropySampling(10)) #


# res = ChainedLinear(neattabX, neattabY)
# bsres = bootstrap(x -> coef(ChainedLinear(x...)), res.data, BasicSampling(50))
# bsres = bootstrap(x -> coef(BraunHolland(x...)), res.data, BasicSampling(50))
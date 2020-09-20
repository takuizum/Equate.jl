using Equate
using Test
using CSV

@testset "Equate.jl" begin
    # ReadData
    include("readdata.jl")
    # SG design
    include("SGdesign.jl")
    # NEAT design
    include("NEATdesign.jl")
    # Smoothing
    include("Smoothing.jl")
end



# Respect from Distributions.jl
println("Potentially stale exports: ")
display(Test.detect_ambiguities(Equate))
println()
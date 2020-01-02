module Equate

using DataFrames, Statistics, Distributions
using Optim#: optimize, BFGS
using GLM#: glm, predict, @formula, Loglink

export
    # Basic function
    freqtab,

    # Equate functions
    Equipercentile,
    Linear,
    Tucker,
    LevineCongeneric,
    ChainedLinear,
    FrequencyEstimation,
    BraunHolland,
    ChainedEquipercentile,

    # Smoothing functions
    presmoothing,
    KernelSmoothing,

    # Support functions
    round2,
    PRF,
    CDF,
    PFu,
    PFl,
    LogLinearFormula,
    ObservableStats,
    BandwidthPenalty,
    EstBandwidth,

    # struct
    ResultEquipercentile,
    ResultLinear,
    ResultTucker,
    ResultLevineCongeneric,
    ResultChainedLinear,
    ResultFrequencyEstimation,
    ResultBraunHolland,
    ResultChainedEquipercentile

abstract type EquateDesign end
"""
    EG
Equivalent (single) group design.
"""
abstract type EG <: EquateDesign end
"""
    NEAT
Non-Equivalent group with Anchor Test design.
"""
abstract type NEAT <: EquateDesign end
abstract type SGEquateMethod end
abstract type NEATEquateMethod end

include("freqtab.jl")
include("SGdesign.jl")
include("NEATdesign.jl")
include("Smoothing.jl")
#-----------------
end # module

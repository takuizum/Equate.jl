module Equate
# =======

using DataFrames
using Statistics: mean, std, cov, cor
using Distributions: Poisson
using Optim: optimize, BFGS
using GLM: glm, predict, @formula, LogLink

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

    # Utility functions
    coef,
    ExpandTable,

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
If the population of 2 tests is identical, the test design will be called the single group design.
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
include("coef.jl")
include("ExpandTable.jl")

"""
    A Julia package for test equating.

Major features:

* Single Group (SG)* design

- `Linear` provides the linear equating.
- `Equipercentile` provides the equipercentile equating.

*Non Equivalent group with Anchor Test (NEAT) design*

- `Tucker`
- `ChainedLinear`
- `ChainedEquipercentile`
- `BraunHolland`
- `FrequencyEstimation`

and more...

According to von Davier, Holland, and Thayer (2004), the test equating has separate 5 steps.

- First step, *Pre-smoothing*.
- Second step is *Estimating the score probabilities*.
- Third step is *Continuization*.
- Fourth step is *Equating*.
- Fifth, last step is *Calculating the Standard Error of Equating(SEE)*.

"""
Equate
# =======
end # module

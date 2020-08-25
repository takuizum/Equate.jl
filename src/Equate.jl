module Equate
# =======

using DataFrames
import Statistics: mean, std, cov, cor, var
import Distributions: Poisson, Normal, pdf, cdf
import Optim: optimize, BFGS
import GLM: glm, predict, @formula, LogLink, coef
import StatsModels: TableRegressionModel
# import RecipesBase: @recipe, plot
using RecipesBase
import Plots: cgrad, @layout
import StatsBase: coef

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
abstract type EquateMethod end
abstract type SGEquateMethod <: EquateMethod end
abstract type NEATEquateMethod <: EquateMethod end

include("freqtab.jl")
include("SGdesign.jl")
include("NEATdesign.jl")
include("Smoothing.jl")
include("coef.jl")
include("ExpandTable.jl")
include("recipe.jl")

"""
    A Julia package for test equating.

Major features:

*Single Group (SG)* design

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
#-----------------
end # module

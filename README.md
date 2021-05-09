# Equate

![CI](https://github.com/takuizum/Equate.jl/workflows/CI/badge.svg)
[![Codecov](https://codecov.io/gh/takuizum/Equate.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/takuizum/Equate.jl)
[![Coveralls](https://coveralls.io/repos/github/takuizum/Equate.jl/badge.svg?branch=master)](https://coveralls.io/github/takuizum/Equate.jl?branch=master)

Equate test scores under the equivalent or non-equivalent group with anchor test design.

## Supported methods

- SG (Single Group design)
  - Linear
  - Mean
  - Equipercentile

- NEAT (Non-Equivalent group design with Anchor Test design)
  - Tucker (Linear equating under some assumptions)
  - Braun & Holland (Linear equating using synthetic population)
  - Chained Linear (also Mean)
  - Chained Equipercentile
  - Frequency Estimation (Equipercentile equating using synthetic population)

- Presmoothing
  - Log linear smoothing with an arbitrary degree.

- Kernel smoothong
  - Gaussian kernel is only supported now.
  - The optimal bandwidth can be estimated.


# Example1: SG design.

1. Prepare data set. Integer or Float vector.
```
using Distributions, Random
Random.seed!(1234)
X = rand(BetaBinomial(100, 4, 10), 500);
Y = rand(BetaBinomial(100, 6, 10), 500);
```

2. Convert the data vector to `FreqTab`

```
# `data` must be Real vector
julia> ftX = freqtab(X; scale = 0:1:100)
Frequency table stats.
         N :      500 
   Missing :        0 
       min :        1 
      maxs :       71 
         μ : 28.39000 
         σ : 12.68921 
  kurtosis : 0.00127 
  skewness : 0.52615 


julia> ftY = freqtab(Y; scale = 0:1:100)
Frequency table stats.
         N :      500 
   Missing :        0 
       min :        6 
      maxs :       77 
         μ : 37.04000 
         σ : 13.11115 
  kurtosis : -0.44152 
  skewness : 0.22442 
```

3. Presmoothing by using `presmoothing`
4. (Optional) Continuization by using `KernelSmoothing`
5. Equate score X to scale Y by the arbitrary method.
```
# Linear Equating
julia> eq_lin = Linear(ftX, ftY)
Equating design: EG
Equated method: Linear.
To show the table, extract `table` element.

# Equipercentile equating
julia> eq_eqp = Equipercentile(ftX, ftY)
Equating design: EG
Equated method: Equipercentile(lower).
To show the table, extract `table` element.
```
6. Evaluate SEE (Standard Error of Equating). Now, Only `BasicSampling(n)` is supported.

```
julia> using Bootstrap, Random
julia> Random.seed!(1234)

julia> @time bootse_lin = bootstrap(x -> coef(Linear(x...)), eq_lin.data, BasicSampling(1000))
  1.727280 seconds (10.12 M allocations: 458.347 MiB, 5.56% gc time, 63.60% compilation time)
Bootstrap Sampling
  Estimates:
     Var │ Estimate  Bias         StdError
         │ Float64   Float64      Float64
    ─────┼─────────────────────────────────
       1 │  1.03325  0.000560431  0.043795
       2 │  7.70597  0.00756638   1.2771
  Sampling: BasicSampling
  Samples:  1000
  Data:     NamedTuple{(:X, :Y), Tuple{Equate.FreqTab, Equate.FreqTab}}: { X 500 × Y 500 }
```


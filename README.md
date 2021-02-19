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


# How to use (SG design)

1. Prepare data set. Integer or Float vector.
2. Convert the data vector to `FreqTab`

```
# `data` is numeric vector
ftX = freqtab(dataX)
dfY = freqtab(dataY)
```

3. Presmoothing by using `presmoothing`
4. (Optional) Continuization by using `KernelSmoothing`
5. Equate score X to scale Y by arbitrary method.
```
Linear(ftX, ftY)
```
6. (Coming soon...) Evaluate SEE.


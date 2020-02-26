# Equate

[![Build Status](https://travis-ci.com/takuizum/Equate.jl.svg?branch=master)](https://travis-ci.com/takuizum/Equate.jl)
[![Codecov](https://codecov.io/gh/takuizum/Equate.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/takuizum/Equate.jl)
[![Coveralls](https://coveralls.io/repos/github/takuizum/Equate.jl/badge.svg?branch=master)](https://coveralls.io/github/takuizum/Equate.jl?branch=master)

Equate test scores under the equivalent or non-equivalent group with anchor test design.

## Supported methods

- SG (Single Group design)
  - Linear
  - Equipercentile

- NEAT (Non-Equivalent group design with Anchor Test design)
  - Tucker (Linear equating under some assumptions)
  - Braun & Holland (Linear equating using synthetic population)
  - Chained Linear
  - Chained Equipercentile
  - Frequency Estimation (Equipercentile equating using synthetic population)

# How to use (SG design)

1. Prepare data set. Integer or Float vector.
2. Convert the data vector to `FreqTab`

```
# `data` is numeric vector
ftX = freqtab(dataX)
dfY = freqtab(dataY)
```
3. Equate score X to scale Y.
4.
```
Linear(ftX, ftY)
```


<!---
# MEMO
using PkgTemplates
tmp = Template(;
    user = "takuizum",
    license = "MIT",
    authors = "Takumi Shibuya",
    dir = "C:\\Users\\bc0089985\\Documents\\GitHub\\",
    julia_version = v"1.3.0",
    ssh = false,
    plugins=[
        TravisCI(),
        Codecov(),
        Coveralls(),
            ]
    )

generate("Equate",tmp)

--->


# Version Update History

### 0.1.2

*A new feature*
- `Linear` and `BraunHolland` functions were updated to return not only the concordance table but also equating coefficients.

- Change specification of arguments of `presmoothing`, `fml`, as the non named arg.

### 0.1.1

Reverse equating direction to match the result to R's `equate` package.

### 0.1.0

A first release.

using Equate

Equate.DiscreteDistribution()


## Haebara
# Symmetric mode
M = Equate.HB(Equate.par_L, Equate.par_K)
M.Forward
Equate._equate!(M)
@code_warntype Equate._equate!(M)
M.slope
M.intercept

# Forward mode
M = Equate.HB(Equate.par_L, Equate.par_K)
M.Reverse.eval = false
Equate._equate!(M)
M.slope
M.intercept

# Reverse mode
M = Equate.HB(Equate.par_L, Equate.par_K)
M.Forward.eval = false
Equate._equate!(M)
M.slope
M.intercept


## SL
# Symmetric mode
M = Equate.SL(Equate.par_L, Equate.par_K)
Equate._equate!(M)
@code_warntype Equate._equate!(M)
M.slope
M.intercept

# Forward mode
M = Equate.SL(Equate.par_L, Equate.par_K)
M.Reverse.eval = false
Equate._equate!(M)
M.slope
M.intercept

# Reverse mode
M = Equate.SL(Equate.par_L, Equate.par_K)
M.Forward.eval = false
Equate._equate!(M)
M.slope
M.intercept
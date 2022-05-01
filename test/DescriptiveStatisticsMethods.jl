using Equate

M = Equate.MM(Equate.par_L, Equate.par_K)
Equate._equate!(M)
M
M.slope
M.intercept

M = Equate.MS(Equate.par_L, Equate.par_K)
Equate._equate!(M)
M
M.slope
M.intercept
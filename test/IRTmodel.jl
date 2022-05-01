using Equate, Random, Distributions

## Blank model

## Logistic model
test = rand(Categorical([0.5, 0.5]), 100_000) .- 1
M = Logistic1(1, 0, 0)
IRF(M, 0.0, 1)
ERF(M, 0.0)
loglik(M, 0.0, test)
@code_warntype loglik(M, 0.0, test)
using Plots
plot(x -> IRF(M, x, 0))
plot(x -> IRF(M, x, 1))


## Graded Model
test = rand(Categorical([0.3, 0.5, 0.2]), 100_000) .- 1
M = Graded(1.0, [1.0, 2.0])
M.b
IRF(M, 0.0, 0)
IRF(M, 0.0, 1)
IRF(M, 0.0, 2)
ERF(M, 0.0)
loglik(M, 0.0, test)
@code_warntype loglik(M, 0.0, test)

using Plots
plot(x -> IRF(M, x, 0))
plot!(x -> IRF(M, x, 1))
plot!(x -> IRF(M, x, 2))

lt!(M, 1.2, 3.0)
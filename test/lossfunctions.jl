using Equate, Distributions

M1 = Graded(1.0, [1.0, 2.0])
M2 = Equate.lt(M1, 1.2, 0.5)

[Equate.HΔ(i, M1, M1) for i in -4:1:4]
[Equate.HΔ(i, M1, M2) for i in -4:1:4]

# using Plots
# plot(x -> IRF(M1, x, 0))
# plot!(x -> IRF(M1, x, 1))
# plot!(x -> IRF(M1, x, 2))
# plot!(x -> IRF(M2, x, 0))
# plot!(x -> IRF(M2, x, 1))
# plot!(x -> IRF(M2, x, 2))

@time Hdiff(0.0, Equate.par_L, Equate.par_K)
@time SLdiff(0.0, Equate.par_L, Equate.par_K)

@code_warntype Hdiff(0.0, Equate.par_L, Equate.par_K)
@code_warntype SLdiff(0.0, Equate.par_L, Equate.par_K)

theta = range(-4, 4; length = 121)
W, Θ = Equate.discrete_normal(theta, 0, 1)
@time Equate.loss1(W, Θ, Equate.par_L, Equate.par_K, 1.2, 0.5)
@code_warntype Equate.loss1(W, Θ, Equate.par_L, Equate.par_K, 1.1, 0.5)

using Optim
# SL
opt = optimize(x -> Equate.SLcrit(W, Θ, Equate.par_L, Equate.par_K, x[1], x[2]), [1.0, 0.0], Newton())
fieldnames(typeof(opt))
opt.minimizer
# HB
opt = optimize(x -> Equate.Hcrit(W, Θ, Equate.par_L, Equate.par_K, x[1], x[2]), [1.0, 0.0], Newton())
fieldnames(typeof(opt))
opt.minimizer
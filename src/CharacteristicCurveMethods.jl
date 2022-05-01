abstract type WeightDistributions end

mutable struct DiscreteDistribution{P <: Real, V <: AbstractVector{P}} <: WeightDistributions
    Θ::V
    W::V
    eval::Bool
end

function DiscreteDistribution(; node = range(-4, 4, length = 121), weight = discrete_normal(node, 0, 1))
    w, n = collect.(weight)
    any(w .< 0) && throw(DomainError("Weights vector must be positive vector."))
    DiscreteDistribution(n, w, true)
end


## Haebara
mutable struct HB{P <: IRTmodels, V <: AbstractVector{P}, R <: Real, D <: WeightDistributions} <: IRTEquateMethod
    To::V
    From::V
    slope::R
    intercept::R
    method::String
    Equated::Bool
    Forward::D
    Reverse::D
    fit
end

"""
    HB(To::V, From::V; Forward::D =DiscreteDistribution(), Reverse::D =DiscreteDistribution()) where {P <: IRTmodels, V <: AbstractVector{P}, D <: WeightDistributions}
Struct for Haebara equating method.
"""
function HB(To::V, From::V; Forward::D =DiscreteDistribution(), Reverse::D =DiscreteDistribution()) where {P <: IRTmodels, V <: AbstractVector{P}, D <: WeightDistributions}
    return HB(To, From, Inf, Inf, "Haebara", false, Forward, Reverse, nothing)
end

function equate!(M::HB; optimize_method = Newton())
    if M.Forward.eval && M.Reverse.eval
        opt = optimize(
            x -> 0.5 * (Equate.Hcrit(M.Forward.W, M.Forward.Θ, M.To, M.From, x[1], x[2]) + Equate.Hcrit(M.Reverse.W, M.Reverse.Θ, M.From, M.To, x[1], x[2]; f = Equate.rlt)), 
            [1.0, 0.0], optimize_method
        )
    elseif M.Forward.eval
        opt = optimize(
            x -> Equate.Hcrit(M.Forward.W, M.Forward.Θ, M.To, M.From, x[1], x[2]), 
            [1.0, 0.0], optimize_method
        )
    elseif M.Reverse.eval
        opt = optimize(
            x -> Equate.Hcrit(M.Reverse.W, M.Reverse.Θ, M.From, M.To, x[1], x[2]; f = Equate.rlt), 
            [1.0, 0.0], optimize_method
        )
    end
    M.slope = opt.minimizer[1]
    M.intercept = opt.minimizer[2]
    M.Equated = true
    M.fit = opt
    return nothing
end

## Sticking-Lord
mutable struct SL{P <: IRTmodels, V <: AbstractVector{P}, R <: Real, D <: WeightDistributions} <: IRTEquateMethod
    To::V
    From::V
    slope::R
    intercept::R
    method::String
    Equated::Bool
    Forward::D
    Reverse::D
    fit
end

function SL(To::V, From::V; Forward::D =DiscreteDistribution(), Reverse::D =DiscreteDistribution()) where {P <: IRTmodels, V <: AbstractVector{P}, D <: WeightDistributions}
    return SL(To, From, Inf, Inf, "Stocking&Lord", false, Forward, Reverse, nothing)
end

function equate!(M::SL; optimize_method = Newton())
    if M.Forward.eval && M.Reverse.eval
        opt = optimize(
            x -> .5 * (Equate.SLcrit(M.Forward.W, M.Forward.Θ, M.To, M.From, x[1], x[2]) + Equate.Hcrit(M.Reverse.W, M.Reverse.Θ, M.From, M.To, x[1], x[2]; f = Equate.rlt)), 
            [1.0, 0.0], optimize_method
        )
    elseif M.Forward.eval
        opt = optimize(
            x -> Equate.SLcrit(M.Forward.W, M.Forward.Θ, M.To, M.From, x[1], x[2]), 
            [1.0, 0.0], optimize_method
        )
    elseif M.Reverse.eval
        opt = optimize(
            x -> Equate.SLcrit(M.Reverse.W, M.Reverse.Θ, M.From, M.To, x[1], x[2]; f = Equate.rlt), 
            [1.0, 0.0], optimize_method
        )
    end
    M.slope = opt.minimizer[1]
    M.intercept = opt.minimizer[2]
    M.Equated = true
    M.fit = opt
    return nothing
end
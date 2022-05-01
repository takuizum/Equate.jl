

## Diff function
function HΔ(θ::L, P1::T, P2::T) where {L <: Real, T <: DichotomousResponseModels}
    SqE(IRF(P1, θ, 1), IRF(P2, θ, 1))
end

function HΔ(θ::L, P1::T, P2::T) where {L <: Real, T <: PolytomousResponseModels}
    P1.K !== P1.K && throw(ArgumentError("The number of category differs between P1 and P2."))
    d = 0.0
    for k in 0:P1.K-1
        d += SqE(IRF(P1, θ, k), IRF(P2, θ, k))
    end
    return d
end

function Hdiff(θ::L, V1::AbstractArray{T, 1}, V2::AbstractArray{T, 1}) where {L <: Real, T <: DiscreteResponseModels}
    length(V1) !== length(V2) && throw(ArgumentError("The length of items differs between V1 and V2."))
    d = 0.0
    for (P1, P2) in zip(V1, V2)
        d += HΔ(θ, P1, P2)
    end
    return d / length(V1)
end

function SLdiff(θ::L, V1::AbstractArray{T, 1}, V2::AbstractArray{T, 1}) where {L <: Real, T <: DiscreteResponseModels}
    length(V1) !== length(V2) && throw(ArgumentError("The length of items differs between V1 and V2."))
    d1, d2 = 0.0, 0.0
    for (P1, P2) in zip(V1, V2)
        d1 += ERF(P1, θ)
        d2 += ERF(P2, θ)
    end
    return SqE(d1 / length(V1), d2 / length(V1))
end

## Loss function

"""
    SLcrit(W::AbstractVector{L}, Θ::AbstractVector{L}, V1::AbstractArray{T, 1}, V2::AbstractArray{T, 1}, A::S, K::S; f::Function = lt) where {L <: Real, T <: DiscreteResponseModels, S <: Real}
    Hcrit(W::AbstractVector{L}, Θ::AbstractVector{L}, V1::AbstractArray{T, 1}, V2::AbstractArray{T, 1}, A::S, K::S; f::Function = lt) where {L <: Real, T <: DiscreteResponseModels, S <: Real}
Objective function for characteristic curve equating methods (Stocking & Lord and Haebara method).
These functions sum `SLdiff` or `Hdiff` along `Θ` vector with the weights `W`. `V1` is a base scale that is target scale of equating, and `V2` is a scale to be converted.
To conduct reverse direction linking, set the base scale as `V2` and specify the named argument `f = rlt`. rlt means 'Reverse Linear Transformation'.


"""
function SLcrit(W::AbstractVector{L}, Θ::AbstractVector{L}, V1::AbstractArray{T, 1}, V2::AbstractArray{T, 1}, A::S, K::S; f::Function = lt) where {L <: Real, T <: DiscreteResponseModels, S <: Real}
    d = 0.0
    V2′ = f.(V2, A, K)
    for (w, θ) in zip(W, Θ)
        d += SLdiff(θ, V1, V2′) * w
    end
    return d
end

function Hcrit(W::AbstractVector{L}, Θ::AbstractVector{L}, V1::AbstractArray{T, 1}, V2::AbstractArray{T, 1}, A::S, K::S; f::Function = lt) where {L <: Real, T <: DiscreteResponseModels, S <: Real}
    d = 0.0
    V2′ = f.(V2, A, K)
    for (w, θ) in zip(W, Θ)
        d += Hdiff(θ, V1, V2′) * w
    end
    return d
end

## Support functions
SqE(x::T, y::T) where {T <: Real} = (x - y)^2

function discrete_normal(M, μ, σ)
    p = pdf.(Normal(μ, σ), M)
    return p ./ sum(p), M
end
abstract type IRTmodels end

abstract type DiscreteResponseModels <: IRTmodels end
abstract type DichotomousResponseModels <: DiscreteResponseModels end
abstract type PolytomousResponseModels <: DiscreteResponseModels end
abstract type OederedPolytomousResponseModels <: PolytomousResponseModels end
abstract type ContinuousResponseModels <: IRTmodels end

Base.copy(x::T) where T <: IRTmodels = T([getfield(x, k) for k in fieldnames(T)]...)

## Elements to define a response model
# 1. Model type, Discrete or Continuous.
# 2. Model (mutable) struct contains item parameters.
# 3. IRF (Item Response Function) that takes an item response as input and, returns a probability as an output.


## Model mutable struct

"""
    Logistic1(a, b, c)
Logistic item response model for dichotomous item response. `Logistic1` stands for the 'discrimination-difficulty' parametrizaion model.

"""
mutable struct Logistic1{A <: Real, B <: Real, C <: Real, R <: Real, M <: AbstractMatrix{R}} <: DichotomousResponseModels
    a::A
    b::B
    c::C
    σ::M
end

"""
    Logistic2(a, b, c)
Logistic item response model for dichotomous item response. `Logistic2` stands for the 'slope-intercept' parametrizaion model.

"""
mutable struct Logistic2{A <: Real, D <: Real, C <: Real, R <: Real, M <: AbstractMatrix{R}} <: DichotomousResponseModels
    a::A
    d::D
    c::C
    σ::M
end

"""
    Graded(a::A, b::AbstractVector{A}) where {A <: Real}
Graded response model (Samejima, 1969), extended dichotomouls IRT model to polytomous one, models the ordered categorical response.
This struct contains parameters for a single item. To represent the whole parameters of a test, put `Graded() into Vector.

- `a` is a discrimination parameter scalar,usually positive.
- `b` is vector of difficulty parameters that sorted in ascending order.
- `σ` is variance-covariance matrix of estimated item parameters.
"""
mutable struct Graded{A <: Real, B <: AbstractVector{A}, R <: Real, M <: AbstractMatrix{R}} <: OederedPolytomousResponseModels
    a::A
    b::B
    K::Int64
    σ::M
end

function Graded(a::A, b::AbstractVector{A}) where {A <: Real}
    !isascending(b) && throw(ArgumentError("Difficulty vector is unordered!"))
    Graded(a, b, length(b) + 1, zeros(length(b) + 1, length(b) + 1))
end

## IRF
function IRF(M::T, θ::Real, x::Int64) where {T <: DichotomousResponseModels}
    s = M.a * (θ - M.b)
    if x == 1
        logistic(s)
    elseif x == 0
        1.0 - logistic(s)
    else
        return 0.0
    end
end

function IRF(M::Graded, θ::Real, x::Int64)
    if x == (M.K - 1)
        p = logistic(M.a * (θ - M.b[end]))
    elseif x == 0
        p = 1.0 - logistic(M.a * (θ - M.b[begin]))
    else
        p = logistic(M.a * (θ - M.b[x])) - logistic(M.a * (θ - M.b[x+1]))
    end
    !(1.0 ≥ p ≥ 0.0) && throw(DomainError("$(p)\nProbability will be (0, 1]. Make sure your parameters or response"))
    return p
end

## Expected IRF
function ERF(M::T, θ::Real) where {T <: DichotomousResponseModels}
    IRF(M, θ, 1)
end

function ERF(M::T, θ::Real) where {T <: OederedPolytomousResponseModels}
    s = 0.0
    for k in 1:M.K-1 # Skip 0 because p×0 is always 0.
        s += k * IRF(M, θ, k)
    end
    return s
end

## Log Likelihood function
function loglik(M::T, θ::Real, x::AbstractVector{Int64}) where {T <: IRTmodels}
    l = 0.0
    for i in x
        l += log(IRF(M, θ, i))
    end
    return l
end

## Linear transformation
function lt(M::Logistic1, A::T, K::T) where {T <: Real}
    N = copy(M)
    N.a = N.a / A
    N.b = A * N.b + K
    return N
end

function lt!(M::Logistic1, A::T, K::T) where {T <: Real}
    M.a = M.a / A
    M.b = A * M.b + K
    return nothing
end

function lt(M::Graded, A::T, K::T) where {T <: Real}
    N = copy(M)
    N.a = N.a / A
    N.b = @. A * N.b + K
    return N
end

function lt!(M::Graded, A::T, K::T) where {T <: Real}
    M.a = M.a / A
    M.b = @. A * M.b + K
    return nothing
end

## Reverse linear transformation
function rlt(M::Logistic1, A::T, K::T) where {T <: Real}
    N = copy(M)
    N.a = M.a * A
    N.b = (M.b - K) / A
    return N
end

function rlt!(M::Logistic1, A::T, K::T) where {T <: Real}
    M.a = M.a * A
    M.b = (M.b - K) / A
    return nothing
end

function rlt(M::Graded, A::T, K::T) where {T <: Real}
    N = copy(M)
    N.a = M.a * A
    N.b = @. (M.b - K) / A
    return N
end

function rlt!(M::Graded, A::T, K::T) where {T <: Real}
    M.a = M.a * A
    M.b = @. (M.b - K) / A
    return nothing
end

## Wrapper

equate()

## Support functions

"""
    isascending(v)
Check if the vector is sorted by ascending.

# Examples
```julia
julia> isascending([0,1,2,3,4])
true

julia> isascending([-1.0, 0.0, 1.0])
true

julia> isascending(reverse([-1.0, 0.0, 1.0]))
false

julia> isascending([-1.0, 0.0, 1.0, 0.0])
false

```

"""
function isascending(v)
    ord = sortperm(v)
    return all(diff(ord) .== 1)
end
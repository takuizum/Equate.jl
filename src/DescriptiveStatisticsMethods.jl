## Mean & Mean methods
mutable struct MM{P <: IRTmodels, V <: AbstractVector{P}, R <: Real} <: IRTEquateMethod
    To::V
    From::V
    slope::R
    intercept::R
    method::String
    Equated::Bool
end

function MM(To::V, From::V) where {P <: IRTmodels, V <: AbstractVector{P}}
    return MM(To, From, Inf, Inf, "Mean & Mean", false)
end

function equate!(M::MM)
    a_to, b_to = extract_a(M.To), extract_b(M.To)
    a_from, b_from = extract_a(M.From), extract_b(M.From)
    # Equate
    A = mean(a_from) / mean(a_to)
    K = mean(b_to) - A * mean(b_from)
    M.slope = A
    M.intercept = K
    M.Equated = true
    return nothing
end

## Mean & Sigma methods
mutable struct MS{P <: IRTmodels, V <: AbstractVector{P}, R <: Real} <: IRTEquateMethod
    To::V
    From::V
    slope::R
    intercept::R
    method::String
    Equated::Bool
end

function MS(To::V, From::V) where {P <: IRTmodels, V <: AbstractVector{P}}
    return MS(To, From, Inf, Inf, "Mean & Sigma", false)
end

function equate!(M::MS)
    b_to = extract_b(M.To)
    b_from = extract_b(M.From)
    # Equate
    A = std(b_to) / std(b_from)
    K = mean(b_to) - A * mean(b_from)
    M.slope = A
    M.intercept = K
    M.Equated = true
    return nothing
end

## Mean & Robust mean methods


## Mean & Geometric mean methods

## Support functions
extract_a(M::T) where {T <: AbstractArray{Logistic1, 1}} = [getfield(m, :a) for m in M]
extract_b(M::T) where {T <: AbstractArray{Logistic1, 1}} = [getfield(m, :b) for m in M]

extract_a(M::T) where {P <: Graded, T <: AbstractArray{P, 1}} = [getfield(m, :a) for m in M]
function extract_b(M::T) where {P <: Graded, T <: AbstractArray{P, 1}}
    v = [getfield(m, :b) for m in M]
    return vcat(v...)
end


function Base.show(O::IO, M::T) where {T <: IRTEquateMethod}
    if M.Equated
        println(O, 
            "$(M.method) Equating\n* slope     = $(M.slope)\n* intercept = $(M.intercept)"
        )
    else
        println(O, 
            "Initialized * $(M.method) Equating"
        )
    end
end
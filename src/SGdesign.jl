# Natural Round
round2(x; digits = 0) = sign(x) * floor( abs(x) * 10.0^digits + 0.5 ) / (10.0^digits)
# Percentile Rank Function
function CDF(x, F::EG)
    if x < minimum(F.table.scale) return 0 end
    if x > maximum(F.table.scale) return 1 end
    F.table.cumprob[F.table.scale .== x][1]
end
function PRF(x, F::EG)
    if x < (minimum(F.table.scale) - F.interval/2.0) return 0.0 end
    if x ≥ (maximum(F.table.scale) + F.interval/2.0) return 100.0 end
    x′ = round2(x)
    Fx1 = CDF(x′-F.interval, F)#F.table.cumfreq[F.table.scale .== (x⃰-1.0)]
    Fx = CDF(x′, F)#F.table.cumfreq[F.table.scale .== x⃰]
    P = 100*(Fx1+(x-x′+F.interval/2.0)*(Fx-Fx1))[1]
    return P
end
# Percentile Function
function p_search_descend(P, F::EG, offset)
    x = nothing;iter = length(F.table.scale)
    while x === nothing
        iter -= 1
        x =  100CDF(F.table.scale[iter], F) > P ? nothing : F.table.scale[iter+offset]
    end
    return x
end
function p_search_ascend(P, F::EG, offset)
    x = nothing;iter = 0
    while x === nothing
        iter += 1
        x = 100CDF(F.table.scale[iter], F) < P ? nothing : iter == 1 ? 0.0 : F.table.scale[iter+offset]
    end
    return x
end
function PFu(P, F::EG)
    if P ≥ 100.0 return (maximum(F.table.scale) + .5) end
    xu = P > 50.0 ? p_search_descend(P, F, 1) : p_search_ascend(P, F, 0)
    x = (P/100 - CDF(xu-F.interval, F)) / (CDF(xu, F) - CDF(xu-F.interval, F))
    return isinf(x) || isnan(x) ? xu -F.interval/2.0 : x + xu -F.interval/2.0
end
function PFl(P, F::EG)
    if P ≤ 0.0 return -.5 end
    xl = P > 50.0 ? p_search_descend(P, F, 0) : p_search_ascend(P, F, -1)
    x = (P/100 - CDF(xl, F)) / (CDF(xl+F.interval, F) - CDF(xl, F))
    return isinf(x) || isnan(x) ? xl + F.interval/2.0 : x + xl + F.interval/2.0
end
# equipercentile equating
struct SGEquateResult <: SGEquateMethod
    method
    table
    estimates
    data
end
"""
    Equipercentile(X::EG, Y::EG; case = :middle)

Equipercentile Equating under equivalent (single) group design. The smoothed frequency can be used.
A table returned by the function contains 2 columns. The first one stands for the base scale score in the test X.
The second one is equated score from test Y, which is equipercentile score on the base scale.

`case` represents which equating case use. Pass to the symbols below.

- `:upper` Calculating scores correspond to arbitrary percentile P, use the smallest integer score with a cumulative percent that is greater than P.
- `:lower` (default) Contrary to above case, use the largest integer score with a cumulative percent that is less than P
- `:both` (Not for practice) Show both case of upper and lower.
- `:middle` Use midpoint case upper between lower.

In my example cases, it seems that equate package in R uses `lower` case to equate.
"""
function Equipercentile(X::EG, Y::EG; case = :lower)
    scaleX = X.table.scale
    eYxu = zeros(Float64, length(scaleX))
    eYxl = zeros(Float64, length(scaleX))
    for (i,v) in enumerate(scaleX)
        P = PRF(v, X)
        eYxu[i] = PFu(P, Y)
        eYxl[i] = PFl(P, Y)
    end
    if case == :upper
        eYx = eYxu
    elseif case == :lower
        eYx = eYxl
    elseif case == :both
        eYx = string.(eYxu, "_", eYxl)
    elseif case == :middle
        eYx = (eYxu .+ eYxl) ./ 2.0
    end
    tbl = DataFrame(scaleX = scaleX, eYx = eYx)
    return SGEquateResult(
        Symbol("Equipercentile($(case))"),
        tbl, 
        nothing, 
        (X = X, Y = Y)
    )
end
# linear equating
"""
    Linear(X::EG, Y::EG)

Linear equating under the equivalent group desing.
This method, equate to match first 2 moments, is so simple to comprehend the equating result.
"""
function Linear(X::EG, Y::EG)
    μX = mean(X.raw); σX = std(X.raw)
    μY = mean(Y.raw); σY = std(Y.raw)
    slope = σY/σX; intercept = μY - slope*μX
    lYx = @. X.table.scale * slope + intercept
    tbl = DataFrame(scaleX = X.table.scale, lYx = lYx)
    return SGEquateResult(
        :Linear,
        tbl, 
        (slope = slope, intercept = intercept), 
        (X = X, Y = Y)
    )
end

# mean equating
"""
    Mean(X::EG, Y::EG)

Mean equating under the equivalent group desing.
This method, equate to match only first moments = mean, is so simple to comprehend the equating result.
"""
function Mean(X::EG, Y::EG)
    μX = mean(X.raw)
    μY = mean(Y.raw)
    slope = 1.0; intercept = μY - μX
    lYx = @. X.table.scale * slope + intercept
    tbl = DataFrame(scaleX = X.table.scale, lYx = lYx)
    return SGEquateResult(
        :Mean, 
        tbl, 
        (slope = slope, intercept = intercept), 
        (X = X, Y = Y)
    )
end

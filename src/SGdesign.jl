# Natural Round
round2(x; digits = 0) = sign(x) * floor( abs(x) * 10.0^digits + 0.5 ) / (10.0^digits)
# Percentile Rank Function
function CDF(x, F::EG)
    if x < minimum(F.tab.scale) return 0 end
    if x > maximum(F.tab.scale) return 1 end
    F.tab.cumprob[F.tab.scale .== x][1]
end
function PRF(x, F::EG)
    if x < (minimum(F.tab.scale) - F.interval/2.0) return 0.0 end
    if x ≥ (maximum(F.tab.scale) + F.interval/2.0) return 100.0 end
    x′ = round2(x)
    Fx1 = CDF(x′-F.interval, F)#F.tab.cumfreq[F.tab.scale .== (x⃰-1.0)]
    Fx = CDF(x′, F)#F.tab.cumfreq[F.tab.scale .== x⃰]
    P = 100*(Fx1+(x-x′+F.interval/2.0)*(Fx-Fx1))[1]
    return P
end
# Percentile Function
function p_search_descend(P, F::EG, offset)
    x = nothing;iter = length(F.tab.scale)
    while x == nothing
        iter -= 1
        x =  100CDF(F.tab.scale[iter], F) > P ? nothing : F.tab.scale[iter+offset]
    end
    return x
end
function p_search_ascend(P, F::EG, offset)
    x = nothing;iter = 0
    while x == nothing
        iter += 1
        x = 100CDF(F.tab.scale[iter], F) < P ? nothing : iter == 1 ? 0.0 : F.tab.scale[iter+offset]
    end
    return x
end
function PFu(P, F::EG)
    if P ≥ 100.0 return (maximum(F.tab.scale) + .5) end
    xu = P > 50.0 ? p_search_descend(P, F, 1) : p_search_ascend(P, F, 0)
    x = (P/100 - CDF(xu-F.interval, F)) / (CDF(xu, F) - CDF(xu-F.interval, F))
    return isinf(x) || isnan(x) ? xu -F.interval/2.0 : x + xu -F.interval/2.0
end
function PFl(P, F::EG)
    if P ≤ 0.0 return -.5 end
    xl = P > 50.0 ? p_search_descend(P, F, 0) : p_search_ascend(P, F, -1)
    x = (P/100 - CDF(xl, F)) / (CDF(xl+F.interval, F) - CDF(xl, F))
    return isinf(x) ? xl + F.interval/2.0 : x + xl + F.interval/2.0
end
# equipercentile equating
struct ResultEquipercentile <: SGEquateMethod
    table::DataFrame
end
"""
    Equipercentile(X::EG, Y::EG; case = :middle)
Equipercentile Equating under equivalent (single) group design. The smoothed frequency can be used.

`case` represents which equating case use. Pass to the symbols below.

- `:upper` Calculating scores correspond to arbitrary percentile P, use the smallest integer score with a cumulative percent that is greater than P.
- `:lower` Contrary to above case, use the largest integer score with a cumulative percent that is less than P
- `:both` (Not for practice) Show both case of upper and lower.
- `:middle` (default) Use midpoint case upper between lower.
"""
function Equipercentile(X::EG, Y::EG; case = :middle)
    scaleY = Y.tab.scale
    eYxu = zeros(Float64, length(scaleY)); eYxl = zeros(Float64, length(scaleY))
    for (i,v) in enumerate(scaleY)
        P = PRF(v, Y)
        eYxu[i] = PFu(P, X)
        eYxl[i] = PFl(P, X)
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
    tbl = DataFrame(scaleY = scaleY, eYx = eYx)
    return ResultEquipercentile(tbl)
end
# linear equating
struct ResultLinear <: SGEquateMethod
    table::DataFrame
end
"""
    Linear(X::EG, Y::EG)
Linear equating under the equivalent group desing.
This method, equate to match first 2 moments, is so simple to comprehension equating result.
"""
function Linear(X::EG, Y::EG)
    μX = mean(X.raw); σX = std(X.raw)
    μY = mean(Y.raw); σY = std(Y.raw)
    slope = σY/σX; intercept = μY - slope*μX
    eYx = @. X.tab.scale * slope + intercept
    tbl = DataFrame(scaleX = X.tab.scale, eYx = eYx)
    ResultLinear(tbl)
end

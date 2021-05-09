
"""
    coef(x::EquateMethod)

Return equating coefficient. The input object is equated result under the SG and NEAT design method.
"""
function coef(x::EquateMethod)
    if isnothing(x.estimates)
        x.table[:, 2]
    else
        x.estimates
    end
end



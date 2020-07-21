
"""
    coef(x::EquateMethod)

Return equating coefficient. The input object is equated result under the SG and NEAT design method.
"""
function coef(x::EquateMethod)
    if all(isnothing.(match.(r"estimates", string.(fieldnames(typeof(x))))))
        println("$(typeof(x)) has no field `estimates`")
    else
        x.estimates
    end
end



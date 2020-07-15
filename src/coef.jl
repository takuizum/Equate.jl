
"""
    coef(x::NEATEquateMethod)

Return equating coefficient. The input object is equated result under the NEAT design method.
"""
function coef(x::NEATEquateMethod)
    if all(isnothing.(match.(r"estimates", string.(fieldnames(typeof(x))))))
        println("$(typeof(x)) has no field `estimates`")
    else
        x.estimates
    end
end


"""
    coef(x::SGEquateMethod)

Return equating coefficient. The input onject is equated result under the SG design method.
"""
function coef(x::SGEquateMethod)
    #
end
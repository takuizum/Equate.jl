using Test, Equate
include("readdata.jl")




function autoscaling(x; float_unit = 0.01)
    if eltype(x) <: Union{Missing, Int64}
        return minimum(x):1:maximum(x)
    elseif eltype(x) <: Union{Real, Missing}
        # get unique elements witout missing values and count them
        returm minimum(x):float_unit:maximum(x)
    elseif eltype{x} <: Union{Missing, String}
        throw(ArgumentError("String is not permitted as variables for freqtab."))
    else
        throw(ArgumentError("Unexpected eltype was passed to `autoscaling`. Please inspect data or set `scale` automatically."))
    end
end

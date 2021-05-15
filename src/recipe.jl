# Check score distributions

function shorten_stats(x)
    txt = @sprintf "N:%i, Mis:%i, Min:%i, Max:%i, μ:%2.2f, σ:%2.2f" x.N x.Missing x.min x.max x.μ x.σ
    return txt
end

@recipe function f(x::FreqTab)
    seriestype := :bar
    linecolor --> :black
    fillcolor --> :transparent
    label --> ""
    title --> shorten_stats(x.stats)
    x.table.scale, x.table.freq
end

@recipe function f(x::SmoothedFreqTab)
    label --> "Smoothed degree = $(length(coef(x.fit))-1)"
    x.table.scale, x.table.freq
end

@recipe function f(x::KernelFreqTab)
    label --> "Kernel smoothed observed probability"
    x.table.scale, x.table.freq
end

# @recipe function f(x::SGFreqTab)
#     # layout
#     layout --> @layout [
#         tophist           _
#         scatter{0.5w,0.5h} righthist
#     ]
#     link := :y
#     framestyle := [:none :axes :axes]
#     legend --> false
#     # scatter plot
#     @series begin
#         seriestype := :scatter
#         markeralpha --> 0.4
#         markershape --> :+
#         markercolor --> :black
#         subplot := 2
#         x.rawX, x.rawV
#     end
#     # Hist : xaxis
#     @series begin
#         seriestype := :line
#         seriescolor := :black
#         subplot := 1
#         x.tableX.scale, x.tableX.freq
#     end
#     # Hist : yaxis
#     @series begin
#         seriestype := :line
#         seriescolor := :black
#         orientation := :h
#         subplot := 3
#         x.tableV.scale, x.tableV.freq
#     end
# end

# Show result of equating
@recipe function f(x::SGEquateResult)
    label --> "Linear"
    xguide --> "Test X"
    yguide --> "Test Y (scaled)"
    legend --> :topleft
    x.table.lYx, x.table.scaleX
end
# NEAT design
@recipe function f(x::NEATEquateResult)
    label --> "Tucker"
    xguide --> "Test X"
    yguide --> "Test Y (scaled)"
    legend --> :topleft
    x.table.lYx, x.table.scaleX
end


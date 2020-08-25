# Check score distributions
@recipe function f(x::FreqTab)
    seriestype := :line
    linecolor --> :gray
    @series begin
        label --> ""
        seriestype := :scatter
        markershape --> :x
        markercolor --> :black
        x.tab.scale, x.tab.freq
    end
    label --> "Empirical probability"
    x.tab.scale, x.tab.freq
end

@recipe function f(x::SmoothedFreqTab)
    label --> "Smoothed degree = $(length(coef(x.fit))-1)"
    x.tab.scale, x.tab.freq
end

@recipe function f(x::KernelFreqTab)
    label --> "Kernel smoothed observed probability"
    x.tab.scale, x.tab.freq
end

# SG design
@recipe function f(x::SGFreqTab)
    label --> "Marginal Histogram under the SG design"
    seriestype := :heatmap
    seriescolor --> cgrad([:white, :gray, :black])
    x.marginal
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
#         x.tabX.scale, x.tabX.freq
#     end
#     # Hist : yaxis
#     @series begin
#         seriestype := :line
#         seriescolor := :black
#         orientation := :h
#         subplot := 3
#         x.tabV.scale, x.tabV.freq
#     end
# end

# Show result of equating
@recipe function f(x::ResultLinear)
    label --> "Linear"
    xguide --> "Test X"
    yguide --> "Test Y (scaled)"
    legend --> :topleft
    x.table.lYx, x.table.scaleX
end

@recipe function f(x::ResultEquipercentile)
    label --> "Equipercentile"
    xguide --> "Test X"
    yguide --> "Test Y (scaled)"
    legend --> :topleft
    x.table.eYx, x.table.scaleX
end

# NEAT design
@recipe function f(x::ResultFrequencyEstimation)
    label --> "Frequency Estimation"
    xguide --> "Test X"
    yguide --> "Test Y (scaled)"
    legend --> :topleft
    x.table.eYx, x.table.scaleX
end

@recipe function f(x::ResultChainedLinear)
    label --> "Chained Linear"
    xguide --> "Test X"
    yguide --> "Test Y (scaled)"
    legend --> :topleft
    x.table.lYx, x.table.scaleX
end

@recipe function f(x::ResultTucker)
    label --> "Tucker"
    xguide --> "Test X"
    yguide --> "Test Y (scaled)"
    legend --> :topleft
    x.table.lYx, x.table.scaleX
end

@recipe function f(x::ResultBraunHolland)
    label --> "Braun & Holland"
    xguide --> "Test X"
    yguide --> "Test Y (scaled)"
    legend --> :topleft
    x.table.lYx, x.table.scaleX
end

@recipe function f(x::ResultChainedEquipercentile)
    label --> "Chained Equipercentile"
    xguide --> "Test X"
    yguide --> "Test Y (scaled)"
    legend --> :topleft
    x.table.eYx, x.table.scaleX
end


@recipe function f(x::SmoothedFreqTab, fit::TableRegressionModel)
    label --> "Smoothed degree = $(length(coef(fit))-1)"
    x.tab.scale, predict(fit)
end

@recipe function f(x::SmoothedFreqTab)
    label --> "Smoothed observed probability"
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
    color --> cgrad([:white, :gray, :black])
    x.marginal
end

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


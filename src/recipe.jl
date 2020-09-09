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
    # seriescolor --> cgrad([:white, :gray, :black])
    x.marginal
end

@recipe function f(x::NEATEquateResult)
    label --> string(x.method)
    xguide --> "Test X"
    yguide --> "Test Y (scaled)"
    legend --> :topleft
    x.table[!, r"Yx"][!, 1], x.table.scaleX
end


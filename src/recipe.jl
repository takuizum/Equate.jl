@recipe function plot(x::SmoothedFreqTab, fit::TableRegressionModel)
    label --> "Smoothed degree = $(length(GLM.coef(fit))-1)"
    (x.tab.scale, predict(fit))
end

@recipe function plot(x::SmoothedFreqTab)
    label --> "Observed probability"
    (x.tab.scale, x.tab.freq)
end
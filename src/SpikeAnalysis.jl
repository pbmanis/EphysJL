module SpikeAnalysis
using LsqFit
using Statistics
using Printf
using FindPeaks1D


function AnalyzeSpikes(tdat, vdat, idat; timewindow=[0.1, 0.6])
    # tdat, vdat, idat are npoints x ntraces
    refractory_ms = 0.001 # msec
    prominence = 0.005
    min_width = 0.00025 # sec
    height = -0.01

    ntraces = size(tdat, 2)
    spikecounts = zeros(Int, ntraces)
    peakamps = Dict{Int, Vector{Float64}}()
    peaklatencies = Dict{Int, Vector{Float64}}()
    for i in 1:ntraces
        dt = mean(diff(tdat[:, 1]))
        pts = findall((tdat[:, i] .>= timewindow[1]) .& (tdat[:, i] .< timewindow[2]))
        vtrace = vdat[pts, i]
        ttrace = tdat[pts, i]
        peaks, properties = findpeaks1d(vtrace; height=height,
            distance=round(Int, refractory_ms / dt),
            width=round(Int, min_width / dt),
            prominence=prominence
        )
        # println("i: ", i, "     peaks: ", peaks)
        spikecounts[i] = length(peaks)
        if length(peaks) > 0
            # println(" vtrace peaks: ", vtrace[peaks])
            # println(" ttrace peaks: ", ttrace[peaks])
            peakamps[i] = vtrace[peaks]
            peaklatencies[i] = ttrace[peaks]
        else
            peakamps[i] = Vector{Float64}()
            peaklatencies[i] = Vector{Float64}()
        end
    end
    return spikecounts, peakamps, peaklatencies
end

end

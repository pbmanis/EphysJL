__precompile__()
module SpikeAnalysis
using LsqFit
using Statistics
using Printf
using FindPeaks1D

#=
This module provides functions for performing analysis of spikes.
1. spike detection: 3 methods:
	a. voltage threshold (with minimum width)
	b. slope and voltage-dependent (Hight and Kalluri)
	c. peak (default)
This module roughly duplicates the ephys.SpikeAnalysis class in python

NOTE: All times are in seconds.
      All voltages are in Volts.
	  All currents are in Amperes.

Data structures are used to simplify handling the results.
Voltagetrace holds the data for a single trace
Spike holds the analysis results for a single spike
Spikes holds an array of individual Spike data
SpikeDetection holds the parameters for the detection algorithm
=#

struct Spike
    indices::Tuple{Int64,Int64}
    trace::Int64 # trace of origin
    eventno::Int64 # spike number in trace
    latency::Float64 # latency in seconds
    amplitude::Float64  # peak kamplitude in mV
    onsettime::Float64
    peaktime::Float64
    risedvdt::Float64
    halfwidth::Float64
    falldvdt::Float64
    ahpdepth::Float64
    mode::String
    detector::String
end

struct Spikes
    label::String
    events::Vector{Spike}
end

mutable struct SpikeDetectParameters
    threshold::Float64 # threshold in V
    refractory::Float64# refractory period in sec
    min_halfwidth::Float64# minimum halfwidth in sec
    HK_dt2::Float64# specific parameters for Hight & Kalluri detector
    HK_C1::Float64 # V
    HK_C2::Float64 # V
    mode::String  # detection mode (schmitt, threshold or peak)
    detector::String  # detector method to use (threshold, argrelmax, Kalluri)
end

#=
aliased function to copy data into the parameters
=#
function SpikeDetectParameters(;
    threshold = 0.0,
    refractory = 0.0007,
    min_halfwidth = 0.0001,
    HK_dt2 = 0.00175,
    HK_C1 = -0.012,
    HK_C2 = 0.011,
    mode = "threshold",
    detector = "Kalluri",
)
    return (SpikeDetectParameters(
        threshold,
        refractory,
        min_halfwidth,
        HK_dt2,
        HK_C1,
        HK_C2,
        mode,
        detector,
    ))
end

"""
Analyze spikes using H&K detection
"""

function AnalyzeSpikes2(
    tdat,
    vdat,
    idat;
    timewindow = [0.1, 0.6],
    pars = Union{SpikeDetectParameters,missing},
)
    # vtrace = VoltageTrace(tdat, vdat, idat, 0)
    # println(pars)
    if ismissing(pars)
        pars = SpikeDetectParameters()
    end
    if pars.detector == "Kalluri"
        spkt, spkv, spki = HightKalluri(tdat, vdat, idat, pars)
        # println(spkt, spki)
    end
    return spkt, spkv, spki
end


function HightKalluri(tdat, vdat, idat, pars)
    #=
    Find spikes using a box method:
    Voltage must be > threshold, and have slope values restricted to a range
    Units must be consistent: x, dt, d2 (s or ms)
    Unist must be consistent: y, thr, C1, C2 (V or mV)
    Note: probably works best with mV and ms, given the constants above.
    to C1, C2 and the width dt2
    From Hight and Kalluri, J Neurophysiol., 2016
    Returns an array of indices in x where spikes occur
    =#
    # parameters for the Hight and Kalluri detecctor
    dt2 = pars.HK_dt2  # s
    C1 = pars.HK_C1  # V
    C2 = pars.HK_C2  # V
    threshold = pars.threshold

    npts = size(vdat)
    spikes = zero(vdat)
    dt = tdat[2] - tdat[1] # sample rate
    iwid = round(Int64, dt2 / dt) # width of region to look for slope values
    for i = iwid:(length(vdat)-iwid)
        if vdat[i] > threshold # when voltage exceeds threshold, check shape
            if (vdat[i] > vdat[i-1]) & (vdat[i] > vdat[i+1])  # local peak
                if ((vdat[i+iwid] - vdat[i]) < C1) & ((vdat[i] - vdat[i-iwid]) > C2)
                    spikes[i] = 1.0
                end
            end
        end
    end
    spki = @view(spikes[2:end]) - @view(spikes[1:end-1]) .> 0
    insert!(spki, 1, 0) # offset for first element
    spike_indices = findall(x -> x == 1, spki)
    spkt = (spike_indices .- 1) .* dt
    spkv = vdat[spike_indices]
    return spkt, spkv, spike_indices
end

"""
AnalyzeSpikes using peak find routine
"""
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

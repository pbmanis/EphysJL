module IVAnalysis
using LsqFit
using Statistics
using Printf
using Revise
using DSP

using PyPlot
using PyCall
pygui(true)
pyimport("matplotlib.pyplot").ion()

include("Acq4Reader.jl")
include("SpikeAnalysis.jl")

export IV_read_and_plot, IV
"""
This module provides functions for computing current-voltage relationships
from step protocols, and for fitting selected traces to a single exponential,
measuring time constants and the sag due to Ih activatoin.

Also included is a simple plotting routine to show the results in a manner
similar to that used in the Python ephys program.

"""


"""
Compute the current-voltage relations of a data set over
both the steady-state window and the peak wincow

tdat, vdat and idate are the time, voltage and current traces 
    (dimensions are npoints, ntraces)
ss_win is the steady-state window
pk_win is the window for measuring the peak voltage deflection
Returns V, I and a minimum value
"""
function compute_iv(tdat, vdat, idat; ss_win=[0.5, 0.6], pk_win=[0.1, 0.2])


    pts = findall((tdat[:, 1] .>= ss_win[1]) .& (tdat[:, 1] .< ss_win[2]))
    pts2 = findall((tdat[:, 1] .>= pk_win[1]) .& (tdat[:, 1] .< pk_win[2]))
    vm = mean(vdat[pts, :], dims=1)
    im = mean(idat[pts, :], dims=1)
    imp = minimum(idat[pts2, :], dims=1)
    return vm, im, imp
end

"""
Fit multiple traces in a current-voltage data set to measure time constants

tdat, vdat and idate are the time, voltage and current traces 
    (dimensions are npoints, ntraces)
ilim : the min and max currents for fitting
p00 : Initial values for the curve fit (DC, amplitude, 1/tau)
iwin : the window over which the data will be fit
"""
function fit_iv(
    tdat::AbstractArray{Float64,2},
    vdat::AbstractArray{Float64,2},
    idat::AbstractArray{Float64,2};
    ilim::AbstractArray{Float64,1}=[-10e-9, -0.05e-9],
    iwin::AbstractArray{Float64,1}=[0.1, 0.2],
    p00::AbstractArray{Float64,1}=[-0.06, -0.01, 200.0],
)

    # pts =  findall((tdat[:,1] .>= window[1]) .& (tdat[:,1] .< window[2]))
    ipts = findall((tdat[:, 1] .>= iwin[1]) .& (tdat[:, 1] .< iwin[2]))
    imn = mean(idat[ipts, :], dims=1)
    imnp = findall((imn .>= ilim[1]) .& (imn .<= ilim[2]))
    nfits = size(imnp)[1]
    npts = size(ipts)[1]
    @printf("# of fits: %d\n", nfits)
    tfit = Array{Float64,2}(undef, npts, nfits)
    vfit = Array{Float64,2}(undef, npts, nfits)
    params = Array{Any,1}(undef, nfits)
    expmodel(t, p) = p[1] .+ p[2] * exp.(-t ./ p[3])
    for i ∈ 1:nfits
        td = tdat[ipts, i]
        p0 = p00
        vd = vdat[ipts, i]
        fit = curve_fit(expmodel, td .- iwin[1], vd, p0)
        params[i] = fit.param
        @printf(
            "Params: DC= %8.2f mV A = %8.2f mV  Tau = %8.3f ms\n",
            params[i][1] * 1e3,
            params[i][2] * 1e3,
            params[i][3] * 1e3
        )
        tfit[:, i] = td
        vfit[:, i] = expmodel(td .- iwin[1], params[i])
    end
    return tfit, vfit, params
end

"""
Do the IV (current clamap current-voltage relationship) analysis.
We present user with a file selector interface; exit if canceled.
Otherwise, we try to do analysis on the selected directory (Acq4 dataset).
"""
function IV()
    filename = ""
    while filename == ""
        filename = open_dialog("Select Dataset Folder", action=GtkFileChooserAction.SELECT_FOLDER)
        println(filename)
        if (filename != "") & isdir(filename)
            IV_read_and_plot(filename=filename, fits=false, ivs=false, analyze_spikes=false)
            filename = ""
        elseif filename == ""
            return
        end
    end
end

function boxcar(x, b::Int=10)
    sizex = size(x)
    osize = [Int64(floor(sizex[1] / b)), sizex[2]]
    odata = Array{Float64,2}(undef, (osize[1], osize[2]))

    for i ∈ 1:size(odata)[2]
        k = 1
        for j ∈ 1:size(odata)[1]
            if k + b > size(x)[1]
                ke = size(x)[1]
            else
                ke = k + b
            end
            odata[j, i] = mean(x[k:ke, i])
            k += b
        end
    end
    return odata
end

"""
Read an HDF5 file, do a little analysis on the data
	-- ivs and fitting
and plot the result
"""
function IV_read_and_plot(; filename::String="", fits::Bool=true, ivs::Bool=true, analyze_spikes::Bool=true,
    decimate::Int=10, maxsize::Int=100000)
    tdat, idat, vdat, data_info = Acq4Reader.read_hdf5(filename)
    print(data_info)
    return
    
    top_lims, bot_lims = Acq4Reader.get_lims(data_info["clampstate"]["mode"])
    @printf(
        "Data lengths: Time=%d  I=%d  V=%d  [# traces = %d]  Duration: %8.3f sec\n",
        size(tdat)[1],
        size(idat)[1],
        size(vdat)[1],
        size(tdat)[2],
        maximum(tdat)[1],
    )
    # if size is > maxsize, decimate by the decimate factor
    if size(tdat)[1] > maxsize
        tdat = boxcar(tdat, decimate)
        idat = boxcar(idat, decimate)
        vdat = boxcar(vdat, decimate)
    end
    # println("Resampled to : ", size(tdat), size(idat), size(vdat)," maxt = ", maximum(tdat[1]), " dt: ", mean(diff(tdat[:, 1], dims=1)))
    if ivs
        vm, imss, imp = IVAnalysis.compute_iv(tdat, vdat, idat)
    end
    if fits
        tfit, vfit, params =
            IVAnalysis.fit_iv(tdat, vdat, idat, iwin=[0.15, 0.2], ilim=[-1e-9, 10e-12])
        tfit2, vfit2, params2 =
            IVAnalysis.fit_iv(tdat, vdat, idat, iwin=[0.2, 0.6], ilim=[-1e-9, 10e-12])

    end
    if analyze_spikes
        spk_counts, spk_amps, spk_latencies = SpikeAnalysis.AnalyzeSpikes(tdat, vdat, idat)
    end

    # now for the plot
    figurename = "IV Analysis"
    fig, axs = plt.subplot_mosaic("""
AC
AD
AE
AF
BG
""", width_ratios=[3, 1], height_ratios=[1, 1, 1, 1, 1], figsize=(8, 10))
    # ion()  # interactive mode on
    # fig, axs = subplots(2, 2, figsize=(8 10))

    fig.suptitle(figurename, fontsize=14, fontweight="bold")


    axl = vec(collect(values(axs)))
    # clean up Plots
    for ax in axl
        ax.grid(false)
        ax.spines["top"].set_visible(false)
        ax.spines["right"].set_visible(false)
        ax.tick_params(direction="out", length=3, width=0.75)
    end
    ntraces = size(tdat)[2]
    voff = collect(range(0, ntraces * 50e-3, ntraces))
    cumulative_offset = 0.0
    for i in 1:ntraces
        
        imean = mean(idat[:, i])
        if imean > 0.0
            cumulative_offset += 50e-3
            vdat[:, i] .+= cumulative_offset
            spk_amps[i] .+= cumulative_offset
        end
    end

    # voltage traces
    axs["A"].plot(
        tdat,
        vdat,
        # xlims = (-inf, tmax),
        # ylims = bot_lims,
        # legend = false,
        linewidth=0.3,
        color="k",
    )
    axs["A"].set_ylabel("Voltage (V)")
    if analyze_spikes
        for i in 1:ntraces
            if spk_counts[i] > 0
                axs["A"].plot(spk_latencies[i], spk_amps[i], "ro", markersize=0.1)
            end
        end
    end
    # curve fits on the voltage traces
    if fits
        axs["A"].plot(tfit, vfit, linewidth=0.5, linestyle="--", color="r")
        axs["A"].plot(tfit2, vfit2, linewidth=0.5, linestyle="--", color="b")
    end
    # current traces
    axs["B"].plot(
        tdat,
        idat,
        # xlims = (-inf, tmax),
        # ylims = top_lims,
        # label = nothing,
        linewidth=0.5,
        color="k",
    )
    # IV (ss and peak)
    if ivs
        axs["D"].plot(
            hcat(vm[1, :], vm[1, :]),
            hcat(imss[1, :], imp[1, :]),
            linestyle="-",
            linewidth=1.0,
            markersize=3,
            marker="o",
        )
    end
    # splitname = splitpath(filename)
    # shortname = joinpath(splitname[(end-3):end]...)

    tight_layout()
    savefig("IV_Analysis.pdf")
    matplotlib.pyplot.close("all")
    print("Show finished")
end

"""
Read an Acq4 HDF5 file, return time, voltage and current traces.
and plot the result to a pdf file.
"""
function test()
    file = "/Volumes/T7_data/NF107Ai32_Het/2022.02.15_000/slice_000/cell_003/CCIV_1nA_max_1s_pulse_000"
    IV_read_and_plot(filename=file, fits=true, ivs=true, analyze_spikes=true)
end

end

__precompile__()
module IVAnalysis
using LsqFit
using Statistics
using Printf
using Revise
using DSP
using Gtk

using PythonCall
ENV["MPLBACKEND"] = "Qt5Agg"
using PythonPlot
pyplot.ion()

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
A data structure to hold the raw traces for the IV data"
mutable because we offset the positions when plotting
"""
mutable struct IVData
	filename::String
	tdat::Array{}
	vdat::Array{}
	idat::Array{}
	data_info::Dict{}
end

# Steady-state and peak measurements
struct IVss
	vm::Array{}
	imss::Array{}
	ivpk::Array{}
end

struct ExponentialFits
    params1::Array{}
	tfit1::Array{}
	vfit1::Array{}
    params2::Array{}
	tfit2::Array{}
	vfit2::Array{}
end

# mutable because we offset the positions when plotting
mutable struct IVSpikes
	counts::Vector{Int64}
	amplitudes::Dict{}
	latencies::Dict{}
    analyze_spikes::Bool
end

"""
Compute the current-voltage relations of a data set over
both the steady-state window and the peak wincow

tdat, vdat and idate are the time, voltage and current traces 
	(dimensions are npoints, ntraces)
ss_win is the steady-state window
pk_win is the window for measuring the peak voltage deflection
Returns V, I and a minimum value
"""
function compute_iv(IV::IVData; ss_win = [0.5, 0.6], pk_win = [0.1, 0.2])

	pts = findall((IV.tdat[:, 1] .>= ss_win[1]) .& (IV.tdat[:, 1] .< ss_win[2]))
	pts2 = findall((IV.tdat[:, 1] .>= pk_win[1]) .& (IV.tdat[:, 1] .< pk_win[2]))
	vm = mean(IV.vdat[pts, :], dims = 1)
	im = mean(IV.idat[pts, :], dims = 1)
	imp = minimum(IV.idat[pts2, :], dims = 1)
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
	IV_Data::IVData;
	ilim::AbstractArray{Float64, 1} = [-10e-9, -0.05e-9],
	iwin::AbstractArray{Float64, 1} = [0.1, 0.2],
	p00::AbstractArray{Float64, 1} = [-0.06, -0.01, 200.0],
)

	# pts =  findall((tdat[:,1] .>= window[1]) .& (tdat[:,1] .< window[2]))

	ipts = findall((IV_Data.tdat[:, 1] .>= iwin[1]) .& (IV_Data.tdat[:, 1] .< iwin[2]))
	imn = mean(IV_Data.idat[ipts, :], dims = 1)
	imnp = findall((imn .>= ilim[1]) .& (imn .<= ilim[2]))

    nfits = size(imnp)[1]
	npts = size(ipts)[1]
	@printf("# of fits: %d\n", nfits)
	tfit = Array{Float64, 2}(undef, npts, nfits)
	vfit = Array{Float64, 2}(undef, npts, nfits)
	params = Array{Any, 1}(undef, nfits)
	expmodel(t, p) = p[1] .+ p[2] * exp.(-t ./ p[3])
	for i ∈ 1:nfits
		td = IV_Data.tdat[ipts, i]
		p0 = p00
		vd = IV_Data.vdat[ipts, i]
		istep = mean(IV_Data.idat[ipts, i]) - mean(IV_Data.idat[1:10, i])
		fit = curve_fit(expmodel, td .- iwin[1], vd, p0;) #  lower=lb, upper=ub)
		params[i] = fit.param
		@printf(
			"Params: DC= %8.2f mV A = %8.2f mV  Tau = %8.3f ms Istep= %8.3f nA\n",
			params[i][1] * 1e3,
			params[i][2] * 1e3,
			params[i][3] * 1e3,
			istep * 1e9
		)
		tfit[:, i] = td
		vfit[:, i] = expmodel(td .- iwin[1], params[i])
	end
	println("returning fit and params")
	return tfit, vfit, imn, params
end

"""
Do the IV (current clamap current-voltage relationship) analysis.
We present user with a file selector interface; exit if canceled.
Otherwise, we try to do analysis on the selected directory (Acq4 dataset).
"""
function IV()
	filename = ""
	while filename == ""
		filename = Gtk.open_dialog("Select Dataset Folder", action = Gtk.GtkFileChooserAction.SELECT_FOLDER)
		println(filename)
		if (filename != "") & isdir(filename)
			IV_read_and_plot(filename = filename, fits = false, ivs = false, analyze_spikes = false)
			filename = ""
		elseif filename == ""
			return
		end
	end
end

function boxcar(x, b::Int = 10)
	sizex = size(x)
	osize = [Int64(floor(sizex[1] / b)), sizex[2]]
	odata = Array{Float64, 2}(undef, (osize[1], osize[2]))

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
function IV_read_and_plot(; filename::String = "", fits::Bool = true, ivs::Bool = true, 
    analyze_spikes::Bool = true,
	decimate::Int = 10, maxsize::Int = 100000)

	tdat, idat, vdat, data_info = Acq4Reader.read_hdf5(filename)
	IV_data = IVData(filename, tdat, vdat, idat,  data_info)

	top_lims, bot_lims = Acq4Reader.get_lims(data_info["clampstate"]["mode"])
	@printf(
		"Data lengths: Time=%d  I=%d  V=%d  [# traces = %d]  Duration: %8.3f sec\n",
		size(IV_data.tdat)[1],
		size(IV_data.idat)[1],
		size(IV_data.vdat)[1],
		size(IV_data.tdat)[2],
		maximum(IV_data.tdat)[1],
	)
	# if size is > maxsize, decimate by the decimate factor
	if size(IV_data.tdat)[1] > maxsize
		IV_data.tdat = boxcar(IV_data.tdat, decimate)
		IV_data.idat = boxcar(IV_data.idat, decimate)
		IV_data.vdat = boxcar(IV_data.vdat, decimate)
	end
	# println("Resampled to : ", size(tdat), size(idat), size(vdat)," maxt = ", maximum(tdat[1]), " dt: ", mean(diff(tdat[:, 1], dims=1)))
	t0 = parse(Float64, IV_data.data_info["MultiClamp1.pulse_params"]["delay"])*1e-3
	t1 = parse(Float64, IV_data.data_info["MultiClamp1.pulse_params"]["duration"])
	println("Pulse start: ", t0, "  duration: ", t1)
	vm, imss, imp = IVAnalysis.compute_iv(
		IV_data,
		ss_win = [t0+0.8*t1, t0+t1], pk_win = [t0+0.001, t0+0.1])
	IV_ss = IVss(vm, imss, imp)
    println("fits: ", fits)
	if fits
		tfit1, vfit1, params1 =
			IVAnalysis.fit_iv(IV_data; iwin = [t0, t0+0.05], ilim = [-1e-9, 10e-12])
		tfit2, vfit2, params2 =
			IVAnalysis.fit_iv(IV_data; iwin = [t0+0.05, t0+t1], ilim = [-1e-9, 10e-12])
		IV_fits = ExponentialFits(
            params1, tfit1, vfit1, params2, tfit2, vfit2)
	else
		IV_fits = nothing
	end
	println("analyze spikes? ", analyze_spikes)
	if analyze_spikes
		counts, amplitudes, latencies = SpikeAnalysis.AnalyzeSpikes(
            IV_data.tdat, IV_data.vdat, IV_data.idat, timewindow = [t0, t0+t1])
		IV_spikes = IVSpikes(counts, amplitudes, latencies, analyze_spikes)
	else
		IV_spikes=nothing
	end

	plot_IV(
		IV_data, IV_fits = IV_fits, IV_spikes = IV_spikes;)
	return
end
function clean_axes(ax)
    ax.grid(false)
    ax.spines["top"].set_visible(false)
    ax.spines["right"].set_visible(false)
    ax.tick_params(direction="out", length=3, width=0.75)
end

function plot_IV(
	IV_data::IVData;
	IV_fits::ExponentialFits,
	IV_spikes::IVSpikes,
)
 
	# now for the plot
	println("Starting figure")
	fig, axs = PythonPlot.pyplot.subplots(2, 2)
	fig.set_size_inches(8, 10)
	fig.suptitle("IV Analysis x", fontsize=14, fontweight="bold")

	# pyplot.tight_layout()
	# pyplot.show()
	# savefig("IV_Analysis.pdf")
	# pyplot.close("all")
	# # println("Figure done")
	# return
	ntraces = size(IV_data.tdat)[2]
	voff = collect(range(0, ntraces * 50e-3, ntraces))
	cumulative_offset = 0.0
	for i in 1:ntraces
		imean = mean(IV_data.idat[:, i])
		if imean > 0.0
			cumulative_offset += 50e-3
			IV_data.vdat[:, i] .+= cumulative_offset
			IV_spikes.amplitudes[i] .+= cumulative_offset
		end
	end

	# voltage traces
	print("plot A")
    clean_axes(axs[1,1])
	axs[1,1].plot(
		IV_data.tdat,
		IV_data.vdat,
		# xlims = (-inf, tmax),
		# ylims = bot_lims,
		# legend = false,
		linewidth = 0.3,
		color = "k",
	)
	axs[1,1].set_ylabel("Voltage (V)")
	if IV_spikes.analyze_spikes
		for i in 1:ntraces
			if IV_spikes.counts[i] > 0
				axs[1,1].plot(IV_spikes.latencies[i], IV_spikes.amplitudes[i], "ro", markersize = 0.1)
			end
		end
	end
	# curve fits on the voltage traces
	if IV_fits != nothing
        print(" plot C")

		# axs[1,1].plot(IV_fits.tfit1, IV_fits.vfit1, "r--", linewidth = 0.5)
		# axs[1,1].plot(IV_fits.tfit2, IV_fits.vfit2, "b--", linewidth = 0.5)
	end
	# current traces
    print(" plot B")
    clean_axes(axs[2,1])
	axs[2,1].plot(
		IV_data.tdat,
		IV_data.idat,
		# xlims = (-inf, tmax),
		# ylims = top_lims,
		# label = nothing,
		linewidth = 0.5,
		color = "k",
	)
	# IV (ss and peak)
	if IV_data != nothing
        print(" plot D")
        clean_axes(axs[2,2])
        axs[2,2].plot(
			hcat(IV_data.vm[1, :], IV_data.vm[1, :]),
			hcat(IV_data.imss[1, :], IV_data.imp[1, :]),
			linestyle = "-",
			linewidth = 1.0,
			markersize = 3,
			marker = "o",
		)
	end
	# splitname = splitpath(filename)
	# shortname = joinpath(splitname[(end-3):end]...)

	tight_layout()
    pyplot.show()
	savefig("IV_Analysis.pdf")
	pyplot.close("all")
	print("Show finished")
end

"""
Read an Acq4 HDF5 file, return time, voltage and current traces.
and plot the result to a pdf file.
"""
function test()
	file = "/Volumes/T7_data/NF107Ai32_Het/2022.02.15_000/slice_000/cell_003/CCIV_1nA_max_1s_pulse_000"
	IV_read_and_plot(filename = file, fits = true, ivs = true, analyze_spikes = true)
end

function test_plot()
	# Create some data
	x = 0:0.1:2π
	y = sin.(x)
	fig, axes = pyplot.subplots(2, 2)
	figurename = "Test Plot"
	fig.suptitle(figurename, fontsize = 16)
	for ax in axes.flat
		# Create a plot using Matplotlib functions
		ax.plot(x, y, label = "sin(x)")
		ax.set_xlabel("x")
		ax.set_ylabel("y")
		ax.set_title("Sine Wave Plot")
		ax.legend()
		ax.grid(true)
	end
	pyplot.tight_layout()
	pyplot.savefig("sine_wave_plot.png")
	# Display the plot (if not in an interactive environment like IJulia)
	pyplot.show()
	println("finis")
end


end

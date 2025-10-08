# __precompile__()
module Acq4Reader
using AddPackage
@add using Statistics
@add using HDF5
# using Plots
# using PyPlot
# using Printf

@add using ThreadsX
@add using ArraysOfArrays
@add using ElasticArrays
@add using Crayons
@add using Distributed
@add using SharedArrays  # important for parallel/looping
# using ProgressMeter
using InteractiveUtils

using PyCall
# pgc = pyimport("pyqtgraph.configfile")  # use python configuration reader...

# Theoretically, the following should work, but we encounter an error:
# "[NameError: name 'array' is not defined]")
# on the scanner target array structure with lots of "array([ 0.00157848,  0.003536  ])"
scriptdir = @__DIR__
pushfirst!(PyVector(pyimport("sys")."path"), scriptdir)
pgc = pyimport("python.configfile")  # use a local python routine to read the config files

export read_hdf5, get_lims, get_stim_times

pyplot()


function get_subdirs(base_dir)
    #
    subdirs = readdir(base_dir)
    print("subdirs in ", base_dir, ": ", join(subdirs, ", "))
    println()
    strings_to_remove = [".DS_Store", ".index", "._.index"]
    filter!(s -> !(s in strings_to_remove), subdirs)
    print("    Found ", length(subdirs), " subdirectories.")
    println("subdirs: ", join(subdirs, ", "))

    return subdirs
end

"""
    read_hdf5(filename)

    Read the protocol data from the specified metaarray file,
    treating it as an HDF5 file (which it is).

    Return the data as 3 2-D arrays (sweep, time) and .
"""
function read_hdf5(filename)
    #=

    New parallel version, 5/24/2021 pbm
    =#

    okfile = isdir(filename)
    if !okfile
        println("file not found:")
        println("    ", filename)
        error()
    end

    # define standard names for devices
    device = "MultiClamp1.ma"
    # device = "Clamp1.ma"
    laser_device = "Laser-Blue-raw.ma"
    photodiode_device = "Photodiode.ma"
    sweeps = get_subdirs(filename)
    idat = Array{Float64}(undef)
    vdat = Array{Float64}(undef)
    tdat = Array{Float64}(undef)
    nsweeps = size(sweeps)[1]
    println("    Nsweeps: ", nsweeps)
    # temporary allocation so we don't lose scope on arrays
    nwave = 1
    s_idat = SharedArray{Float64,2}((nwave, nsweeps))
    s_vdat = SharedArray{Float64,2}((nwave, nsweeps))
    s_tdat = SharedArray{Float64,2}((nwave, nsweeps))

    first = true
    mode = ""
    clampstate = ""
    data_info = ""
    println("    Number of threads: ", Base.Threads.nthreads())
    # note we set this up for threading, but that causes memory errors...
    # kept same structure here though.

    for s = 1:nsweeps
        time, data, data_info = read_one_sweep(filename, sweeps[s], device)
  
        if time == false
            continue
        end
        # sweep mode will be VC, IC or IC=0 (may depend on age of acquisition code)
        sweep_mode = String(data_info["clampstate"]["mode"])
        if first
            mode = sweep_mode
        else
            if mode != sweep_mode
                throw(
                    ErrorException(
                        string(
                            "Mode changed from ",
                            _mode,
                            " to ",
                            sweep_mode,
                            "inside protocol",
                        ),
                    ),
                )
            end
        end

        indx1, indx2 = get_indices(data_info)

        if first  # we don't know allocation size until reading the first run
            nwave = size(data)[1]
            s_idat = SharedArray{Float64,2}((nwave, nsweeps))
            s_vdat = SharedArray{Float64,2}((nwave, nsweeps))
            s_tdat = SharedArray{Float64,2}((nwave, nsweeps))
            first = false
        end
        s_tdat[:, s] = time
        s_idat[:, s] = data[:, indx1]
        s_vdat[:, s] = data[:, indx2]

    end
    if mode == ""
        finalize(s_tdat)
        finalize(s_idat)
        finalize(s_vdat)
        @everywhere GC.gc()
        println("mode: ", mode)
        throw(ErrorException("No acquisition mode found"))
    end
    idat = deepcopy(s_idat)
    tdat = deepcopy(s_tdat)
    vdat = deepcopy(s_vdat)
    indexfile = joinpath(filename, ".index")
    # println("Reading index file: ", indexfile)
    cf = Configfile.read_configfile(indexfile)
    if haskey(cf["."]["devices"], "Laser-Blue-raw")
        wavefunction = cf["."]["devices"]["Laser-Blue-raw"]["channels"]["pCell"]["waveGeneratorWidget"]["function"]
        data_info["Laser.wavefunction"] = deepcopy(wavefunction)
    elseif haskey(cf["."]["devices"], "MultiClamp1")
        wavefunction = cf["."]["devices"]["MultiClamp1"]["waveGeneratorWidget"]["function"]
        data_info["MultiClamp1.wavefunction"] = deepcopy(wavefunction)
        pulse = cf["."]["devices"]["MultiClamp1"]["waveGeneratorWidget"]["stimuli"]["Pulse"]
        data_info["MultiClamp1.pulse_pars"] = [
            parse(
                Float64,
                pulse["start"]["value"],
            ),
            parse(
                Float64,
                pulse["length"]["value"],
            ),
        ]

    else
        wavefunction = nothing
    end
    finalize(s_tdat)
    finalize(s_idat)
    finalize(s_vdat)
    finalize(wavefunction)
    @everywhere GC.gc()
    println("    Finished reading protocol data.")
    return tdat, idat, vdat, data_info
end

function get_lims(mode)
    #=
    Set the top and bottom limits to some defaults according to the 
    acquistion mode
    =#
    if mode == "'VC'"
        # print("VC")
        println(GREEN_FG, "Limts for VC", WHITE_FG)
        top_lims = (-0.4e-9, 0.4e-9)
        bot_lims = (-120e3, 100e3)
    elseif mode == "'IC'"
        println(GREEN_FG, "Limits for IC", WHITE_FG)
        top_lims = (-120e-3, 40e-3)
        bot_lims = (-2e-9, 2e-9)
    else
        println(RED_FG, "Unknown Mode: ", mode, WHITE_FG)
        top_lims = (-inf, +inf)
        bot_lims = (-inf, +inf)
    end
    return top_lims, bot_lims
end

function get_indices(data_info)
    #= 
    get the array indices that correpond to v and i
    depending on the acquisition mode
    =#
    mode = String(data_info["clampstate"]["mode"])  # will be VC, IC or IC=0 (may depend on age of acquisition code)
    c0_units = String(data_info["c0"]["units"])
    c1_units = String(data_info["c1"]["units"])
    if mode == "'VC'"
        # println("VCMode")
        if c0_units == "'A'"
            topidx = 1
            botidx = 2
        elseif c0_units == "'V'"
            topidx = 2
            botidx = 1
        end
    elseif mode == "'IC'"
        # println("ICMode")
        if c0_units == "'V'"
            topidx = 2
            botidx = 1
        elseif c0_units == "'A'"
            topidx = 1
            botidx = 2
        end
    else
        println("mode is not known: ", mode)
    end
    return topidx, botidx
end

function get_stim_arg(name, data_info; device::AbstractString="Laser")
    names = Dict("latency" => 1, "duration" => 2, "amplitude" => 3)
    index = names[name]
    query = device * ".wavefunction"
    wv = data_info[query]
    u = split(wv, "\\n")
    stim_arg = Vector{Float64}()
    re_float = r"([+-]?\d+\.?\d*), ([+-]?\d+\.?\d*), ([+-]?\d+\.?\d*)"
    for i in 1:size(u)[1]
        s = match(re_float, u[i])
        append!(stim_arg, parse(Float64, s[index]))
        end
    return stim_arg
end



function read_one_sweep(filename::AbstractString, sweep_dir, device)
    #=
    Get one waveform sweep from the protocol for
    the given sweep dir (protocol sequence) and the 
    specified device
    Returns the time array, the sweep data, and some info
    =#

    full_filename = joinpath(filename, sweep_dir, device)
    okfile = isfile(full_filename)
    if !okfile
        println("File not found:")
        println("    ", full_filename)
        return false, false, false
    end
    tf = HDF5.ishdf5(full_filename)
    if !tf
        return false, false, false
    end
    fid = h5open(full_filename, "r")

    c0 = h5readattr(full_filename, "info/0/cols/0")
    c1 = h5readattr(full_filename, "info/0/cols/1")
    c2 = h5readattr(full_filename, "info/0/cols/1")
    ClampState = h5readattr(full_filename, "info/2/ClampState")
    DAQPrimary = h5readattr(full_filename, "info/2/DAQ/primary")
    DAQSecondary = h5readattr(full_filename, "info/2/DAQ/secondary")
    DAQCommand = h5readattr(full_filename, "info/2/DAQ/command")

    data_info = Dict(
        "clampstate" => deepcopy(ClampState),
        "c0" => c0,
        "c1" => c1,
        "c2" => c2,
        "DAQ.Primary" => deepcopy(DAQPrimary),
        "DAQ.Secondary" => deepcopy(DAQSecondary),
        "DAQ.Command" => deepcopy(DAQCommand),
        "Laser.wavefunction" => "",
    )

    time_array = deepcopy(fid["info"]["1"]["values"][:])
    # println("    Read sweep: ", sweep_dir, "  time length: ", size(time_array)[1],)
    # println(" dt: ", mean(diff(time_array)))
    data_array = deepcopy(fid["data"][:, :])
    close(fid)

    return time_array, data_array, data_info
end

function test_configread()
    filename = "/Users/pbmanis/Desktop/2018.09.27_000/ImageSequence_000/.index"
    data = pgc.readConfigFile(filename)
    println(data)
    fn = "/Users/pbmanis/Desktop/Python/mrk-nf107/data_for_testing/CCIV/.index"
    data = pgc.readConfigFile(fn)
    println(data)
end
    
function test_reader()
    # local test file name
    # filename = "/Volumes/Pegasus/ManisLab_Data3/Kasten_Michael/Pyramidal/2018.02.12_000/slice_001/cell_000/CCIV_1nA_max_000"
    filename="/Volumes/T7_data/NF107Ai32_Het/2022.03.18_000/slice_000/cell_000/CCIV_4nA_max_1s_pulse_posonly_000"
    # filename = "/T7_data/NF107Ai32_Het/2022.03.18_000/slice_000/cell_000/CCIV_1nA_max_1s_pulse_000"
    read_hdf5(filename)
end

end

function julia_main()::Cint
    # do something based on ARGS?
    return 0 # if things finished successfully
end

module Acq4Reader
using Statistics
using HDF5
using Printf
using Base.Threads
using ArraysOfArrays
using ElasticArrays
using Crayons.Box
using Distributed
using SharedArrays  # important for parallel/looping

include("Configfile.jl")

export read_hdf5, get_lims, get_stim_times
"""
Find the subdirectories of the given base_dir, 
and return those that are protocol sweeps,
excluding .DS_Store and .index files.
"""
function get_subdirs(base_dir)
    # read the subdirectories of the current base_dir
    # only return those that are protocol sweeps
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
    Read the protocol data from the specified metaarray file,
    treating it as an HDF5 file (which it is)

    Return the data, and some "standard" y limits.
    New parallel version, 5/24/2021 pbm
    added parsing of pulse waveforms 10/2025 pbm

"""
function read_hdf5(filename)

    device = "MultiClamp1.ma"
    # device = "Clamp1.ma"
    laser_device = "Laser-Blue-raw.ma"
    photodiode_device = "Photodiode.ma"
    sweeps = get_subdirs(filename)
    idat = Array{Float64}(undef)
    vdat = Array{Float64}(undef)
    tdat = Array{Float64}(undef)
    nsweeps = size(sweeps)[1]
    println(WHITE_FG, "    Nsweeps: ", nsweeps)
    # temporary allocation so we don't lose scope on arrays
    nwave = 1
    s_idat = SharedArray{Float64,2}((nwave, nsweeps))
    s_vdat = SharedArray{Float64,2}((nwave, nsweeps))
    s_tdat = SharedArray{Float64,2}((nwave, nsweeps))

    first = true
    mode = ""
    clampstate = ""
    data_info = ""
    println("    Number of threads: ", Threads.nthreads())
    # note we wet this up for threading, but that causes memory errors...
    # kept same structure here though.

    for s = 1:nsweeps
        time, data, data_info = read_one_sweep(filename, sweeps[s], device)

        if time == false
            continue
        end
        # will be VC, IC or IC=0 (may depend on age of acquisition code)

        sweep_mode = String(data_info["clampstate"]["mode"])
        if first
            mode = sweep_mode
        else
            if mode != sweep_mode
                throw(
                    ErrorException(
                        "Mode changed from ",
                        _mode,
                        " to ",
                        sweep_mode,
                        "inside protocol",
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
            # s_tdat[:,s] = time
            # s_idat[:,s] = data[:, indx1]
            # s_vdat[:,s] = data[:, indx2]
        end
        s_tdat[:, s] = time
        s_idat[:, s] = data[:, indx1]
        s_vdat[:, s] = data[:, indx2]
        # println("s: ", s)
    end
    if mode == ""
        finalize(s_tdat)
        finalize(s_idat)
        finalize(s_vdat)
        @everywhere GC.gc()
        throw(ErrorException("No acquisition mode found"))
    end
    idat = deepcopy(s_idat)
    tdat = deepcopy(s_tdat)
    vdat = deepcopy(s_vdat)
    indexfile = joinpath(filename, ".index")
    cf = Configfile.readConfigFile(indexfile)[2]["."]
    # println(keys(cf["devices"]))
    if haskey(cf["devices"], "Laser-Blue-raw")
        wavefunction =
            cf["devices"]["Laser-Blue-raw"]["channels"]["pCell"]["waveGeneratorWidget"]["function"]
    else
        wavefunction = nothing
    end
    # print("Wavefunction: ", wavefunction, "\n")
    # println(keys(cf["devices"]["MultiClamp1"]["waveGeneratorWidget"]))
    if haskey(cf["devices"]["MultiClamp1"], "waveGeneratorWidget")
        wavefunction_mc =
            cf["devices"]["MultiClamp1"]["waveGeneratorWidget"]["function"]
        # print("MC Wavefunction: ", wavefunction_mc, "\n")
        # parse a string that looks like: pulse(150*ms, 1*s, Pulse_amplitude)
        rp1 = r"(?P<mode>[a-z]+\()(?P<delay>\d+.?\d)\*(?P<delay_unit>[msunp]{0,2})\s*,\s*"
        rp2 = r"(?P<duration>\d*\.?\d*)\*(?P<duration_unit>[msunp]{0,2})\s*,\s*"
        rp3 = r"(?P<Param>[_a-zA-Z]*)\)"
        s = match(rp1*rp2*rp3, wavefunction_mc)
        # println("MC pulse parameters: ", s.captures)
        print("MC Params: ", cf["devices"]["MultiClamp1"]["waveGeneratorWidget"]["params"], "\n")
        data_info["MultiClamp1.pulse_params"] = s
    end
    if haskey(cf["devices"], "Photodiode")
        wavefunction_pd =
            cf["devices"]["Photodiode"]["channels"]["pCell"]["waveGeneratorWidget"]["function"]
        # print("PD Wavefunction: ", wavefunction_pd, "\n")
    end
    if haskey(cf["devices"], "stimuli")
        wavefunction_stimuli =
            cf["devices"]["stimuli"]
        # print("Stimuli Wavefunction: ", wavefunction_stimuli, "\n")
    end
    data_info["Laser.wavefunction"] = deepcopy(wavefunction)
    println("Data info keys: ", keys(data_info) , "\n")
    finalize(s_tdat)
    finalize(s_idat)
    finalize(s_vdat)
    finalize(wavefunction)
    @everywhere GC.gc()
    println("    Finished reading protocol data.")
    return tdat, idat, vdat, data_info
end

    #=
    Set the top and bottom limits to some defaults according to the 
    acquistion mode
    =#
function get_lims(mode)

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

    #= 
    get the array indices that correpond to v and i
    depending on the acquisition mode
    =#
function get_indices(data_info)

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
        println(RED_FG, "mode is not known: ", mode, WHITE_FG)
    end
    return topidx, botidx
end

# function get_stim_times(data_info)
#     wv = data_info["Laser.wavefunction"]
#     u = split(wv, "\n")
#     stim_lats = Vector{Float64}()
#     re_float = r"[+-]?\d+\.?\d*"
#     for i = 1:size(u)[1]
#         s = match(re_float, u[i])
#         append!(stim_lats, parse(Float64, s.match) * 1e3)
#     end
#     return stim_lats
# end

    #=
    Retrieve stimulus times from a device's wavefunction parameters
    The default will get the information from a Laser device
    =#
function get_stim_times(data_info; device="Laser")

    query = device * ".wavefunction"
    wv = data_info[query]
    u = split(wv, "\n")
    stim_lats = Vector{Float64}()
    re_float = r"[+-]?\d+\.?\d*"
    for i = 1:size(u)[1]
        s = match(re_float, u[i])
        append!(stim_lats, parse(Float64, s.match) * 1e3)
    end
    return stim_lats
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
        println(RED_FG, "File not found:")
        println("    ", full_filename, WHITE_FG)
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

"""
Check the configuration file reader by reading a real acq4 index file
"""
function test_configread()
    # filename = "/Users/pbmanis/Desktop/2018.09.27_000/ImageSequence_000/.index"
    # data = Configfile.readConfigFile(filename)
    # println(data)
    # fn = "/Users/pbmanis/Desktop/Python/mrk-nf107/data_for_testing/CCIV/.index"
    # data = Configfile.readConfigFile(fn)
    # println(data)
    file = "/Volumes/T7_data/NF107Ai32_Het/2022.02.15_000/slice_000/cell_003/CCIV_1nA_max_1s_pulse_000"
    # file = "/Volumes/T7_data/NF107Ai32_Het/2022.03.18_000/slice_000/cell_000/CCIV_4nA_max_1s_pulse_posonly_000"
    file1 = file * "/.index"
    data = Configfile.readConfigFile(file1)
    print("\n",file)
    println(data)  
    file2 = file * "/000/.index" 
    data2= Configfile.readConfigFile(file2)  # the sweep index
    print("\n", file2)
    println(data2)

 end

"""
Test the reader by reading a real acq4 data file
""" 
function test_reader()
    # local test file name
    # file= "/Volumes/Pegasus/ManisLab_Data3/Kasten_Michael/Pyramidal/2018.02.12_000/slice_001/cell_000/CCIV_1nA_max_000"
    file = "/Volumes/T7_data/NF107Ai32_Het/2022.03.18_000/slice_000/cell_000/CCIV_4nA_max_1s_pulse_posonly_000"
    # file= "/T7_data/NF107Ai32_Het/2022.03.18_000/slice_000/cell_000/CCIV_1nA_max_1s_pulse_000"
    file = "/Volumes/T7_data/NF107Ai32_Het/2022.02.15_000/slice_000/cell_003/CCIV_1nA_max_1s_pulse_000"
    read_hdf5(file)
end


end

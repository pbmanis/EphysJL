#if isinteractive()
using Pkg
Pkg.activate(".")
# using CondaPkg
# CondaPkg.resolve()
# CondaPkg.add("pyqt")
# Pkg.add("PythonPlot")
# Pkg.instantiate()
# Pkg.precompile()
# include("src/Configfile.jl")
# include("src/Acq4Reader.jl")
# Acq4Reader.test()
# Acq4Reader.test_reader()
include("src/IVAnalysis.jl")
IVAnalysis.test()
# include("src/SpikeAnalysis.jl")

# IVAnalysis.test()
# include("src/MiniAnalysis.jl")
# include("src/LSPSAnalysis.jl")
# include("test/test_pyplot.jl")
#end
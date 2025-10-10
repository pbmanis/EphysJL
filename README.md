# EphysJL

Julia versions of some of the Python ephys analysis routines.
These are generally faster than the python versions, but
also currently very stripped-down

Includes:

* IVAnalysis: Basic IV analysis (current-voltage relationships, firing rates, basic
spike detection, time constants)

* Laser scanning photostimulation analysis, including event detection, plotting traces
 and maps.



Install:
This package needs PyCall. The best way to get this cleanly is to:
remove prior Julia installs, including:
~/.julia
~/.juliup
~/.JuliPro

curl -fsSL https://install.julialang.org | sh 
source ~/.zshrc (or whatever shell you are using)

Go to the directory where your packages are located (mine are at ~/Desktop/Julia)
cd Julia/EphysJL

julia -L Main.jl

Then, in julia,
>julia Acq4Reader.test_reader()
>julia IVAnalysis.test()



enter Julia

* >julia> using Pkg;
* julia> ]
* pkg> activate .  # must be done in the project directory because that is where the env is defined
* julia> include("IVAnalysis.jl") # must do before including your program
* julia> include("src/Acq4Reader.jl")
* julia>Acq4Reader.test_configure()
.. should print out a long dict with the configuration information.
* julia>Acq4Reader.test() 
.. should read a test data set (you need to put an IV from acq4 into the test_data directory first).
* julia>IVAnalysis.test()
.. Runs IV and plots the data using pycall/matplotlib.pyplot, to a file (IV_Analysis.pdf" in the main directory)```


May have to do: Pkg.instantiate() at some point as well.

--------------
Delete and rebuild: If you suspect a corrupted Python environment, deleting the .CondaPkg directory and rebuilding can resolve the issue.
Exit Julia.
Navigate to your project directory.
Run rm -rf .CondaPkg.
Restart Julia and run import Pkg; Pkg.build("PythonCall") and import CondaPkg; CondaPkg.resolve(). 

------------------
using PythonCall
using Pkg
Pkg.add("PythonPlot")
ENV["MPLBACKEND"] = "Qt5Agg" # Set the environment variable for PythonPlot
using PythonPlot


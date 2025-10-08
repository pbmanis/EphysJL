# EphyJL

Julia versions of some of the ephys analysis routines. Faster than the python versions, but
stripped-down also

Install:
This package needs PyCall. The best way to get this cleanly is to:
remove prior Julia installs, including:
~/.julia
~/.juliup
~/.JuliPro

curl -fsSL https://install.julialang.org | sh 
source ~/.zshrc (or whatever shell you are using)

cd Julia/EphysJL

Go to the directory where your packages are located (mine are at ~/Desktop/Julia)
enter Julia
>julia> using Pkg;
>julia> ]
pkg> activate .  # must be done in the project directory because that is where the env is defined
julia> include("Revise") # must do before including your program
julia> include("src/Acq4Reader.jl")

julia>Acq4Reader.test_configure()
.. should print out a long dict with the configuration information.



push!(pyimport("sys")["path"], "EphysJL/python")
configr = pyimport("configfile") # no extension

Pkg.add(path="./EphysJL")

May have to do: Pkg.instantiate()

pkg> develop ../LocalPackage


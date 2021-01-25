# Installation of EFTfitter.jl

EFTfitter.jl is written in the [Julia programming language](https://julialang.org/). 
To use it, you need to install Julia itself, the EFTfitter.jl and BAT.jl Julia packages as well as a few additional Julia packages. 


## Installing Julia

Julia is available for Linux, OS-X and Windows, and easy to install:

* [Download Julia](https://julialang.org/downloads/).

* Extract the archive, resp. run the installer.

* You may want to add the Julia `bin` directory to your `$PATH`. To get the location of the Julia `bin` directory on OS-X or Windows, start a Julia session (via applications menu) and run the Julia command `Sys.BINDIR`.



## Installing EFTfitter.jl
EFTfitter.jl is not yet a registered Julia package but still easy to install.
To install it, simply run:
```julia
julia> using Pkg
julia> pkg"add https://github.com/Cornelius-G/EFTfitter.jl#dev"
```
### Installing dependencies: BAT.jl and further Julia packages 
BAT.jl (v.2.0) is required for EFTfitter. To install it, run:
```julia
julia> using Pkg
julia> pkg"add BAT"
```

You will also need to install the following Julia packages to run the EFTfitter.jl examples:
```julia
julia> using Pkg
julia> pkg"add Distributions IntervalSets Plots ValueShapes StatsBase"
```
## Further information
If you'd like to [precompile](https://docs.julialang.org/en/v1/manual/modules/index.html#Module-initialization-and-precompilation-1) all installed packages right aways (otherwise they'll get precompiled when loaded for the first time), run:

```julia
julia> pkg"precompile"
```

For further infortmation on how to set up Julia and BAT.jl see the [installation guide](https://bat.github.io/BAT.jl/dev/installation/) of the [BAT.jl documentation](https://bat.github.io/BAT.jl/dev/).

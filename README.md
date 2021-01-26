
# EFTfitter.jl

[![Build Status](https://github.com/invenia/PkgTemplates.jl/workflows/CI/badge.svg)](https://github.com/tudo-physik-e4/EFTfitter.jl/actions)
[![Documentation for stable version](https://img.shields.io/badge/docs-stable-blue.svg)](https://tudo-physik-e4.github.io/EFTfitter.jl/stable)
[![Documentation for development version](https://img.shields.io/badge/docs-dev-blue.svg)](https://tudo-physik-e4.github.io/EFTfitter.jl/dev)


New implementation of the [EFTfitter](https://github.com/tudo-physik-e4/EFTfitterRelease) in the [Julia languange](https://julialang.org/). 
Tool for constraining the parameters of physics models using Bayesian inference by combining measurements of (different) observables.
Particularly suited for EFT (effective field theory) interpretations. 

Work-in-progress, interfaces and functionalities might be subject to changes.

## Installation
The EFTfitter.jl package can be installed using:
```julia
julia> using Pkg
julia> pkg"add https://github.com/tudo-physik-e4/EFTfitter.jl#master"
```

Please see the [installation guide](https://tudo-physik-e4.github.io/EFTfitter.jl/dev/installation/) for further information on the installation of Julia, EFTfitter.jl and its dependencies.


## Documentation & Tutorials
Please see the [documentation](https://tudo-physik-e4.github.io/EFTfitter.jl/dev/) of EFTfitter.jl for further information.\
Tutorials and examples on how to use EFTfitter.jl can be found [here](https://github.com/tudo-physik-e4/EFTfitter.jl/tree/main/examples/tutorials).


## Template
A template for starting the implementation of your own analysis is provided [here](https://github.com/tudo-physik-e4/EFTfitter.jl/tree/main/examples/empty_template).


## Citing EFTfitter.jl
When using EFTfitter.jl for your work, please consider citing:

Nuno Castro, Johannes Erdmann, Cornelius Grunwald, Kevin Kroeninger, Nils-Arne Rosien, *EFTfitter - A tool for interpreting measurements in the context of effective field theories*,  [Eur. Phys. J. C 76 (2016) 8, 432](https://link.springer.com/article/10.1140/epjc/s10052-016-4280-9)
```
@article{EFTfitter2016,
    author = {Castro, Nuno and Erdmann, Johannes and Grunwald, Cornelius and Kr\"oninger, Kevin and Rosien, Nils-Arne},
    title = "{EFTfitter---A tool for interpreting measurements in the context of effective field theories}",
    eprint = "1605.05585",
    archivePrefix = "arXiv",
    primaryClass = "hep-ex",
    doi = "10.1140/epjc/s10052-016-4280-9",
    journal = "Eur. Phys. J. C",
    volume = "76",
    number = "8",
    pages = "432",
    year = "2016"
}
```

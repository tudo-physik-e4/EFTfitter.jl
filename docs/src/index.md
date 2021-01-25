# EFTfitter.jl

EFTfitter is a tool for constraining the parameters of physics models using Bayesian inference 
by combining measurements of (different) observables. 
It is particularly suited for EFT (effective field theory) models.\
The Bayesian inference is performed using the [*Bayesian Analysis Toolkit - BAT.jl*](https://github.com/bat/BAT.jl).\
A detailed explanation of the EFTfitter approach can be found in the [EFTfitter paper](https://link.springer.com/article/10.1140/epjc/s10052-016-4280-9). 
    
EFTfitter.jl is a new implementation of the previous [C++ version of EFTfitter](https://github.com/tudo-physik-e4/EFTfitterRelease).

Work in progress: Interfaces and functions might be subject to changes.

----
### Short description
For further information on the statistical approach see the
[EFTfitter paper](https://link.springer.com/article/10.1140/epjc/s10052-016-4280-9").

**Assumption**: Measurements of physical quantities are approximately gaussian.   
This allows to combine the measurements using the following likelihood:

```math
p(\vec{x}|\vec{y}) = \sum_{i=1}^{n}\sum_{j=1}^{n} [\vec{x} - U\vec{y}]_i \mathcal{M}_{ij}^{-1} [\vec{x} - U\vec{y}]_j
```
where ``\vec{x}`` is the vector of measurements with length ``n `` and ``\vec{y}`` is the vector of predicitions for the values of the ``N`` observables. 
The element ``U_{ij}`` of the (n x N)-matrix ``U`` is unity if ``x_i`` is a measurement of the observable ``y_j``, and zero otherwise.
The covariance matrix
```math
\mathcal{M}_{ij} = \text{cov}[x_i, x_j] = \sum_{k=1}^{M} \text{cov}^{(k)}[x_i, x_j]
```
is the sum over the covariances of all ``k`` types of uncertainties.

### Required user inputs
* **Parameters**: Model parameters to be fitted.
* **Observables**: Predictions for physical quantities as a function of the free model parameters.
* **Measurements**: Measured values of the observables and uncertainty values for different types of uncertainties.
* **Correlations**: Correlation matrices for the uncertainties.

## Features of EFTfitter.jl
* combine multiple measurements of the same observable
* combine measurements of different observables
* fit model parameters using Bayesian inference
* include (multiple categorizes of) uncertainties and correlations
* use simple input formats for measurements, uncertainties & correlations
* rank the individual measurements and uncertainty types by their impact on the result of the combination
* treat unknown correlation coefficients as nuisance parameters in the fit

## Documentation & Tutorials
Please see the [Tutorial](https://cornelius-g.github.io/EFTfitter.jl/dev/tutorial/) and
[Advanced Tutorial](https://cornelius-g.github.io/EFTfitter.jl/dev/advanced_tutorial/) for
detailed information on how to use EFTfitter. 
Executable versions of the tutorials can be found [here](https://github.com/Cornelius-G/EFTfitter.jl/tree/dev/examples).
You can also try EFTfitter.jl right now at [binder](https://github.com/Cornelius-G/EFTfitter.jl/tree/dev/examples/notebooks).

All functions and types are documented in the
[API Documentation](https://cornelius-g.github.io/EFTfitter.jl/dev/api/).
---
## Citing EFTfitter

When using EFTfitter.jl for your work, please consider citing:

Nuno Castro, Johannes Erdmann, Cornelius Grunwald, Kevin Kroeninger, Nils-Arne Rosien, EFTfitter - A tool for interpreting measurements in the context of effective field theories, [Eur. Phys. J. C 76 (2016) 8, 432](https://link.springer.com/article/10.1140/epjc/s10052-016-4280-9)

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

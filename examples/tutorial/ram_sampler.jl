import BAT: AbstractSamplingAlgorithm, bat_sample_impl
using BAT: DensitySampleVector, AbstractRNG, AbstractMeasureOrDensity, MCMCSampleID, AnyMeasureOrDensity, ConvergenceTest, transform_and_unshape, bat_convergence
using Parameters
using ValueShapes
using ArraysOfArrays
using InverseFunctions

# simple wrapper for the RAM sampler [1011.4381] using RobustAdaptiveMetropolisSampler.jl
using RobustAdaptiveMetropolisSampler

export RAMSampler

@with_kw struct RAMSampler{
    CT<:ConvergenceTest
} <: AbstractSamplingAlgorithm
    trafo = PriorToGaussian()
    nchains::Int = 4
    nsteps::Int = 10^5
    nburnin::Int = floor(Int,0.1*nsteps) # number of burnin steps to throw away per chain
    convergence::CT = BrooksGelmanConvergence()
    x0 = fill(nothing, nchains)
    M0 = fill(1., nchains)
    opt_α = fill(0.234, nchains)
    γ=fill(2/3, nchains)
    q=fill(Normal(), nchains)
end


function bat_sample_impl(
    rng::AbstractRNG,
    target::AnyMeasureOrDensity,
    algorithm::RAMSampler
)
    density_notrafo = convert(BAT.AbstractMeasureOrDensity, target)
    density, trafo = transform_and_unshape(algorithm.trafo, density_notrafo)
    shape = varshape(density)

    ram_vector = []
    @info "Start RAM sampling of $(algorithm.nchains) chains with $(Threads.nthreads()) threads."

    Threads.@threads for i in 1:algorithm.nchains
        try
            push!(ram_vector, generate_ram_samples(i, rng, density, algorithm; show_progress = ifelse(i==1, true, false)))
        catch err
            @info "RAM sampling failed for thread $(Threads.threadid())."
        end
    end

    samples, M, acceptance_rate = collect_samples(ram_vector)

    samples_trafo = shape.(reduce(vcat, samples))
    samples_notrafo = inverse(trafo).(samples_trafo)

    #converged = bat_convergence(samples_notrafo, algorithm.convergence).converged
    #@info "RAM chains have$(ifelse(converged, "", " not")) converged"

    return (result = samples_notrafo, result_trafo = samples_trafo, trafo = trafo, M = M, acceptance_rate = acceptance_rate)

end


function generate_ram_samples(
    i::Integer,
    rng::AbstractRNG,
    target_density::Any,
    algorithm::RAMSampler;
    show_progress = false
)
    n_tries = 0
    success = false

    ram_result_vec = []
    while !success && n_tries < 4000
        n_tries += 1
        try
            x0 = any(isnothing.(algorithm.x0)) ? bat_initval(rng, target_density, InitFromTarget()).result : algorithm.x0[i]
            M0 = algorithm.M0[i]
            println("n_tries: ", n_tries)

            ram_result = RobustAdaptiveMetropolisSampler.RAM_sample(logdensityof(target_density), x0, M0, algorithm.nsteps;
                opt_α=algorithm.opt_α[i], γ=algorithm.γ[i], q=algorithm.q[i], show_progress=show_progress, output_log_probability_x=true)

            success = true
            push!(ram_result_vec, ram_result)
        catch e
           success = false
        end
    end

    ram_result = ram_result_vec[1]
    nburnin = algorithm.nburnin
    samples = nestedview(ram_result.chain')[nburnin+1:end]
    logvals = ram_result.log_probabilities_x[nburnin+1:end]

    return (samples=samples, logvals=logvals, M=ram_result.M, acceptance_rate=ram_result.acceptance_rate)
end


function collect_samples(ram_vector)
    M = getproperty.(ram_vector, :M)
    acceptance_rate = getproperty.(ram_vector, :acceptance_rate)
    samples_raw = getproperty.(ram_vector, :samples)
    logvals_raw = getproperty.(ram_vector, :logvals)
    samples = []

    for i in 1:length(ram_vector)
        n_samples = length(samples_raw[i])
        sample_id = fill(MCMCSampleID(i, 0, 0, 0), n_samples)
        push!(samples, DensitySampleVector(samples_raw[i], logvals_raw[i], info=sample_id))
    end

    return (samples=samples, M=M, acceptance_rate=acceptance_rate)
end

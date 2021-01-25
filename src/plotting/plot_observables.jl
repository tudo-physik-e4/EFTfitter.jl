# # plot recipe for plotting an Observable object
@recipe function f(
    observable::Observable,
    params::Any;
    titlewidth=50
)
    x_nt, x_range, x_i = range_to_namedtuples(params)
    x_label = keys(params)[x_i]
    title_string = get_obs_title(params, x_i, titlewidth=titlewidth)
    y = observable.func.(x_nt)

    y_label = string(observable.func)

    @series begin
        seriestype --> :path
        label --> y_label*"("*string(x_label)*")"
        xguide --> string(x_label)
        yguide --> y_label
        title --> title_string
        x_range, y
    end

end


@recipe function f(
    observable::Observable,
    param::NamedTuple,
    default_params::Any;
    titlewidth=50
)
    params = merge(default_params, param)

    @series begin
        titlewidth --> titlewidth
        observable, params
    end
end

# plot recipe for plotting a particular Observable from a NamedTuple
@recipe function f(
    observables::NamedTuple{<:Any, <:Tuple{Vararg{Observable}}},
    name::Symbol,
    params::Any;
    titlewidth=50
)
    @series begin
        titlewidth --> titlewidth
        observables[name], params
    end
end

# plot recipe for plotting a particular Observable from a NamedTuple
@recipe function f(
    observables::NamedTuple{<:Any, <:Tuple{Vararg{Observable}}},
    name::Symbol,
    params::Any,
    default_params::Any;
    titlewidth=50
)
    @series begin
        titlewidth --> titlewidth
        observables[name], params, default_params
    end
end

# plot recipe for plotting an Array{Observable}
@recipe function f(
    observables::AbstractVector{Observable},
    params::Any;
    titlewidth=50,
    bin_names = ["bin$i" for i in 1:length(observables)]
)
    nbins = length(observables)
    ys = [obs.func(params) for obs in observables]

    @series begin
        seriestype --> :steppre
        label --> "prediction"
        xticks --> (collect(1:nbins).+0.5, bin_names)
        linewidth --> 1.5
        1:nbins+1, [ys[1], ys...]
    end
end

# internal: convert a NamedTuple that contains a range into a vector of NamedTuples
function range_to_namedtuples(params)
    ks = keys(params)
    vals = collect(values(params))

    i = 0
    for j in 1:length(ks)
        if !isa(vals[j], Real)
            i = j
        end
    end

    xs = collect(vals[i])
    nt = Vector{NamedTuple}(undef, length(xs))

    for j in 1:length(xs)
        vals[i] = xs[j]
        nt[j] = namedtuple(collect(ks), vals)
    end

    return nt, xs, i
end


function get_obs_title(params, i; titlewidth=50)
    ks = keys(params)
    vals = collect(values(params))

    title_string = ""
    for j in 1:length(ks)
        if j != i
            title_string *= string(ks[j], "=", vals[j], ", ")
        end
    end
    title_string = title_string[1:end-2]

    return replace(title_string, Regex("(.{$(titlewidth)})") => Base.SubstitutionString("\\1\n"))
end

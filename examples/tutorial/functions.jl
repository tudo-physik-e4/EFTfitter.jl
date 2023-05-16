using BenchmarkTools

struct Prediction
    pred::Float64
    unc::Float64
end

Prediction(a::Float64) = Prediction(a, 0.)
Prediction(a::Tuple{Float64, Float64}) = Prediction(a[1], a[2])
Prediction(a::Prediction) = a

observable1(params) = rand()
observable2(params) = (rand(), rand())

functions_vector_1 = fill(observable1, 100)
functions_vector_2 = fill(observable2, 100)
functions_vector_12 = vcat(fill(observable1, 50), fill(observable2, 50))


# using tuple:
to_tuple(f) = x -> (f(x), 0.)
functions_vector_12_tuple = vcat(fill(to_tuple(observable1), 50), fill(observable2, 50))


# using Prediction struct:
observable1_p(params) = Prediction(rand())
observable2_p(params) = Prediction(rand(), rand())

functions_vector_2_p = fill(observable2_p, 100)
functions_vector_12_p = vcat(fill(observable1_p, 50), fill(observable2_p, 50))


function eval_functions_a(f_vec, params)
    res = [f(params) for f in f_vec]
    return res
end

function eval_functions_a2(f_vec, params)
    res::Vector{Tuple{Float64, Float64}} = [f(params) for f in f_vec]
    return res
end


function eval_functions_c(f_vec, params, preds, uncs)
    for i in eachindex(f_vec)
        res::Tuple{Float64, Float64} = f_vec[i](params)
        preds[i] = res[1]
        uncs[i] = res[2] 
    end
    return preds, uncs
end

function eval_functions_d(f_vec, params, preds, uncs)
    for i in eachindex(f_vec)
        res::Prediction = f_vec[i](params)
        preds[i] = res.pred
        uncs[i] = res.unc
    end
    return preds, uncs
end

function eval_functions_e(f_vec, params, preds, uncs)
    for i in eachindex(f_vec)
        res::Prediction = Prediction(f_vec[i](params))
        preds[i] = res.pred
        uncs[i] = res.unc
    end
    return preds, uncs
end

params = rand(3)

# A: Current EFTfitter Implementation
@btime eval_functions_a(functions_vector_1, params)  
@btime eval_functions_a(functions_vector_2, params)  
@btime eval_functions_a(functions_vector_12, params)
@btime eval_functions_a(functions_vector_12_tuple, params) 

@code_llvm eval_functions_a(functions_vector_12, params) 

# B: 


# C: Using preallocation
preds = zeros(100)
uncs = zeros(100)
@btime eval_functions_c(functions_vector_2, params, preds, uncs) 
@btime eval_functions_c(functions_vector_12_tuple, params, preds, uncs) 


# C: Using preallocation & Prediction Type
preds = zeros(100)
uncs = zeros(100)
@btime eval_functions_d(functions_vector_2_p, params, preds, uncs) 
@btime eval_functions_d(functions_vector_12_p, params, preds, uncs) 

@btime eval_functions_e(functions_vector_12, params, preds, uncs)


@code_warntype eval_functions_d(functions_vector_12_p, params, preds, uncs)

typeof((1.,2.))
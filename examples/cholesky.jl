using LinearAlgebra

function generate_covariance_matrix(size::Int)
    # Generate a random square matrix of size x size
    A = rand(size, size)

    # Construct the covariance matrix by multiplying A with its transpose
    covariance_matrix = A * A'

    return covariance_matrix
end




function cholesky_product(M, x)
    L = cholesky(M).L
    z = L' * x
    return dot(z, z)
end

function cholesky_product2(L, x)
    z = L * x
    return dot(z, z)
end

function my_product(M, x)
    z = M * x
    return dot(x, z)

end




using BenchmarkTools


d = 2000

M = generate_covariance_matrix(d)
x = rand(d)

@btime my_product(M, x)

L = cholesky(M).L';
@btime cholesky_product2(L, x)


# For Model Uncertainties
using BenchmarkTools
using LinearAlgebra

using LinearAlgebra

# Define the lower Cholesky factor L and matrix A
L = [3.0 0.0 0.0; 2.0 1.0 0.0; 1.0 2.0 2.0]
A = [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0]

# Compute M = L * L'
M = L * L'

# Solve the linear system L * X = A
X = L \ A

# Calculate the inverse of M + A using the Woodbury matrix identity
inv_MplusA = M \ (A * M \ I)

# Print the result
println(inv_MplusA)




M + A
@btime inv(M + A)


# For Model Uncertainties
using LinearAlgebra

function inverse_sum(M_inv, D)
    n = size(M_inv, 1)  # Assuming M_inv and D are square matrices of the same size
    
    # Calculate the intermediate matrices
    A = inv(M_inv * D + I(n))
    B = M_inv * D * A
    C = M_inv - B * M_inv
    
    return C
end


# Define the matrices M, D, and M_inv
d = 500
L = rand(d)
M = L*L'  # Example symmetric positive definite matrix
D = Diagonal(rand(d))  # Example diagonal matrix
M_inv = inv(M)  # Inverse of M

# Calculate the inverse of M + D using M_inv
@btime inverse_sum(M_inv, D)

@btime inv($(M + D))


using WoodburyMatrices

W = Woodbury(M, I(d), D, I(d))

Wi = inv(W)
Matrix(Wi)
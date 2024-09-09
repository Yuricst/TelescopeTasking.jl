"""Implementation of greedy policy for single telescope tasking problem"""


"""
    solve_greedy(A, T; exclude_target_indices::Union{Vector{Int},Nothing} = nothing)

Get allocations for single telescope tasking problem via greedy policy
"""
function solve_greedy(A, T; exclude_target_indices::Union{Vector{Int},Nothing} = nothing)
    tstart = time()
    # iterate through chronologically ordered passes and choose ones that is the next earliest
    n,m = size(A)
    Y = zeros(Bool, n)
    X = zeros(Bool, m)

    if isnothing(exclude_target_indices)
        exclude_target_indices = Int[]
    end

    # always use first pass
    Y[1] = true
    X[findfirst(x -> x > 0.5, A[1,:])] = true
    for i in 2:n
        i_last = findlast(x -> x == true, Y)
        k = findfirst(x -> x > 0.5, A[i,:])
        if (T[i_last,i] > 0.5) && !(k in exclude_target_indices)
            Y[i] = true
            X[k] = true
            push!(exclude_target_indices, k)
        end
    end
    tend = time()
    telapsed = tend - tstart
    return X, Y, telapsed
end


"""
    solve_greedy(problem::TelescopeTaskingProblem; exclude_target_indices::Union{Vector{Int},Nothing} = nothing)

Get allocations for single telescope tasking problem via greedy policy
"""
function solve_greedy(problem::TelescopeTaskingProblem; exclude_target_indices::Union{Vector{Int},Nothing} = nothing)
    return solve_greedy(problem.A, problem.T; exclude_target_indices = exclude_target_indices)
end


"""
    solve_greedy(problem::MultiTelescopeTaskingProblem)

Get allocations for multi telescope tasking problem via greedy policy
"""
function solve_greedy(problem::MultiTelescopeTaskingProblem)
    tstart = time()

    # initialize storage
    X = zeros(Int, problem.m)
    Y_per_telescope = Vector[]
    exclude_target_indices = Int[]

    # iterate through telescope
    for (A,T) in zip(problem.A_per_telescope, problem.T_per_telescope)
        Xk,Yk,_ = solve_greedy(A, T; exclude_target_indices)
        push!(Y_per_telescope, Yk)
        X += Xk

        # append indices that have already been observed
        exclude_target_indices = findall(x -> x > 0.5, X)
    end
    
    tend = time()
    telapsed = tend - tstart
    Y = vcat(Y_per_telescope...)
    return X, Y, Y_per_telescope, telapsed
end
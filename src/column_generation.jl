"""Column generation scheme for MTTP sharing passes"""


mutable struct ReducedMasterProblem
    P::Int                  # number of schedule patterns
    s::Int                  # number of telescopes
    B::Matrix               # map matrix from targets to patterns, m-by-P
    w::Vector               # priority coefficient on each target
    num_exposure::Int       # number of exposures required for each target

    objective_value::Real   # objective value

    patterns::Vector        # schedule patterns
    z_val::Vector           # value of z variables chosen
end


function Base.show(io::IO, problem::ReducedMasterProblem)
    println(io, "Reduced master problem for MTTP")
    m,P = size(problem.B)
    selected_z = findall(>(0.5), problem.z_val)
    n_passes = [sum(pat) for pat in problem.patterns[selected_z]]
    println(io, "    Number of telescopes s = $(problem.s)")
    println(io, "    Number of targets m    = $(m)")
    println(io, "    Number of patterns P   = $(P)")
    println(io, "    Pattern indices        = $(selected_z)")
    println(io, "    Number of passes       = $(n_passes)")
end


function solve!(problem::ReducedMasterProblem, solver)
    # extract size of matrices
    m,P = size(problem.B)
    # create model
    model = Model(solver; add_bridges = false)
    set_silent(model)
    @variable(model, X[1:m], Bin);
    @variable(model, z[1:P], Bin);
    @constraint(model, sufficient_exposure[k = 1:m],
        sum(problem.B[k,p] * z[p] for p = 1:problem.P) >= problem.num_exposure * X[k]
    )
    @constraint(model, sum(z) == problem.s)
    @objective(model, Max, sum(problem.w[k] * X[k] for k = 1:m))
    optimize!(model)
    
    # store solutions
    @assert is_solved_and_feasible(model)
    problem.z_val = value.(z)
    X_BLP = value.(X)
    problem.objective_value = objective_value(model)

    # relax integrality constraint of pattern vector and solve to get duals
    unset_binary.(X)
    set_lower_bound.(X, 0.0)
    set_upper_bound.(X, 1.0)
    unset_binary.(z)
    set_lower_bound.(z, 0.0)
    set_upper_bound.(z, 1.0)
    optimize!(model)
    println("   Relaxed RMP objective = ", objective_value(model))

    # Obtain a new dual vector
    termination_status(model)
    @assert is_solved_and_feasible(model; dual = true)
    duals = dual.(sufficient_exposure)
    # @show dual.(sufficient_exposure)
    # reduced_costs = [problem.w[k] - (sum(problem.B[k,p] for p = 1:P) - problem.num_exposure)*duals[k] for k = 1:m]
    reduced_costs = [problem.num_exposure*duals[k] - problem.w[k] for k = 1:m]

    # @show value.(X)
    # @show value.(z)
    # # @show reduced_costs
    # @show reduced_cost.(X)
    # @show reduced_cost.(z)
    for k = 1:m
        if value(X_BLP[k]) > 1e-8
            duals[k] = 0.0
        end
    end
    return duals
end


function master_problem_to_schedules(problem::ReducedMasterProblem)
    Y_per_telescope = Vector[]
    for q = 1:problem.s
        if problem.z_val[q] > 0.5
            push!(Y_per_telescope, problem.patterns[q])
        end
    end
    return Y_per_telescope
end


mutable struct PricingProblem
    duals_iter::Vector
    duals::Vector               # dual variables for each target
    A::Matrix                   # binary n-by-m matrix associating observation arcs to targets
    T::Matrix                   # binary n-by-n matrix for transition feasibility from arc i to j

    objective_value::Real       # objective value
end


function Base.show(io::IO, problem::PricingProblem)
    println(io, "Pricing problem for MTTP")
    n,m = size(problem.A)
    println(io, "    Required exposure E       = $(problem.num_exposure)")
    println(io, "    Number of targets m       = $(m)")
    println(io, "    Number of passes n        = $(n)")
    println(io, "    Current pricing objective = $(problem.objective_value)")
end


function solve!(problem::PricingProblem, solver)
    # extract size of matrices
    n,m = size(problem.A)
    # create model & solve
    model = Model(solver; add_bridges = false)
    set_silent(model)
    @variable(model, X[1:m], Bin);
    @variable(model, Y[1:n], Bin);
    @constraint(model, sufficient_observation[k = 1:m],
                sum(problem.A[i,k] * Y[i] for i = 1:n) >= X[k])
    @constraint(model, transition_feasibility[i=1:n-1, j=i+1:n], 
                Y[i] + Y[j] <= 1 + problem.T[i,j])
    @objective(model, Max,
               sum(problem.duals[k] * X[k] for k=1:m))
    optimize!(model)
    problem.objective_value = objective_value(model)
    
    @assert is_solved_and_feasible(model)
    return value.(Y)
end


function schedule_to_B_column(y::Vector, A::Matrix)
    return [sum(A[i,k]*y[i] for i in 1:size(A,1)) for k in 1:size(A,2)]
end


"""
Solve the multi-telescope tasking problem for the case where passes are shared between telescopes.
"""
function solve_column_generation(
    original_problem::MultiTelescopeTaskingProblem,
    solver,
    initial_patterns::Vector;
    maxiter::Int = 2,
)
    # number of telescopes
    s = length(original_problem.n_per_telescope)        # number of telescopes
    for q = 1:s-1
        @assert original_problem.A_per_telescope[q] == original_problem.A_per_telescope[q+1]
    end

    # initialize reduced master problem
    P = length(initial_patterns)
    B = zeros(Int, original_problem.m, P)
    for (p,y) in enumerate(initial_patterns)
        B[:,p] .= schedule_to_B_column(y, original_problem.A_per_telescope[1])
    end

    reduced_master_problem = ReducedMasterProblem(
        P,
        s,
        B,
        original_problem.w,
        original_problem.num_exposure,
        0.0,                                # place-holder for objective value
        initial_patterns,
        zeros(P),                               # place-holder for z
    )

    # initialize pricing problem
    pricing_problem = PricingProblem(
        Vector[],
        zeros(original_problem.m),
        original_problem.A_per_telescope[1],
        original_problem.T_per_telescope[1],
        1.0,                                # place-holder for objective value
    )

    for iter = 1:maxiter
        # solve reduced master problem
        duals = solve!(reduced_master_problem, solver)
        pricing_problem.duals .= duals
        push!(pricing_problem.duals_iter, duals)

        # solve pricing problem
        new_schedule = solve!(pricing_problem, solver)
        
        # store new schedule
        reduced_master_problem.P += 1
        reduced_master_problem.B = hcat(reduced_master_problem.B, schedule_to_B_column(new_schedule, original_problem.A_per_telescope[1]))

        @printf("Column generation iter %d ... RMP objective = %1.2f with %d columns; Pricing objective = %1.2f; %d/%d positive duals\n",
            iter, reduced_master_problem.objective_value, reduced_master_problem.P, pricing_problem.objective_value,
            sum(duals .> 0), length(duals))

        if pricing_problem.objective_value < 1 + 1e-8
            @printf("Pricing problem objective = %1.2f under 1; optimal solution found!\n",
                pricing_problem.objective_value)
            break
        end

        if new_schedule in reduced_master_problem.patterns
            @printf("Pricing problem objective = %1.2f; new schedule already in patterns; breaking!\n",
                pricing_problem.objective_value)
            break
        else
            push!(reduced_master_problem.patterns, new_schedule)
        end
    end
    return reduced_master_problem, pricing_problem
end
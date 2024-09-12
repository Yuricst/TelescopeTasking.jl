"""Lagrangean Relaxation"""



mutable struct LagrangeanRelaxedMultiTelescopeTaskingProblem
    mu::Vector          # Lagrange multipliers
    n_per_telescope::Vector{Int}        # number of observation arcs per telescope
    A_per_telescope     # binary n-by-m matrix associating observation arcs to targets
    T_per_telescope     # binary n-by-n matrix for transition feasibility from arc i to j
    w::Vector                           # priority coefficient on each target per telescope
    num_exposure::Int                   # number of exposures required for each target
    m::Int              # number of targets
    s::Int              # number of telescopes
    Z_LR::Real          # Lagrangean relaxation objective
    X_LR::Vector        # binary m-vector for each telescope
    Y_LR::Vector        # binary n-by-m matrix for each telescope

    Z_LH::Real          # Lagrangean heuristic objective
    X_LH::Vector        # binary m-vector for each telescope
    Y_LH::Vector        # binary n-by-m matrix for each telescope

    Theta::Real         # step size hyperparameter for subgradient method
end


function Base.show(io::IO, problem::LagrangeanRelaxedMultiTelescopeTaskingProblem)
    println(io, "Lagrange-Relaxation of Multi-Telescope Tasking Problem")
    println(io, "    Number of telescopes s = $(problem.s)")
    println(io, "    Number of targets m    = $(problem.m)")
    println(io, "    Upper-bound objective  = $(problem.Z_LR)")
    println(io, "    Lower-bound objective  = $(problem.Z_LH)")
end


function solve_relaxed_problem!(problem::LagrangeanRelaxedMultiTelescopeTaskingProblem, solver)
    
    # reset objective value
    problem.Z_LR = 0.0

    # solve for X's
    for k = 1:problem.m
        if problem.w[k] >= problem.mu[k] * problem.num_exposure
            problem.X_LR[k] = 1
            problem.Z_LR += (problem.w[k] - problem.mu[k] * problem.num_exposure)
        else
            problem.X_LR[k] = 0
        end
    end

    # solve for Y's
    Y_LR = Vector[]
    for (q,(Aq,Tq)) in enumerate(zip(problem.A_per_telescope, problem.T_per_telescope))
        # construct problem
        nq, m = size(Aq)
        model =  Model(solver, add_bridges = false)
        set_silent(model)
        @variable(model, Yq[1:nq], Bin)
        @constraint(model, transition_feasibility[i=1:nq-1, j=i+1:nq], 
            Yq[i] + Yq[j] <= 1 + Tq[i,j])
        @objective(model, Max, sum(problem.mu[k] * sum(Aq[i,k]*Yq[i] for i = 1:nq) for k = 1:m))
        optimize!(model)
        @assert is_solved_and_feasible(model)

        # store
        problem.Z_LR += objective_value(model)
        Y_val = value.(Yq) .> 0.5;
        push!(Y_LR, Y_val)
    end

    # store back into problem struct
    problem.Y_LR = Y_LR
    return
end


function heuristic_feasible!(problem::LagrangeanRelaxedMultiTelescopeTaskingProblem)
    for k = 1:problem.m
        if sum(sum(problem.A_per_telescope[q][i,k] * problem.Y_LR[q][i] 
                   for i = 1:problem.n_per_telescope[q]) for q = 1:problem.s) >= problem.num_exposure
            problem.X_LH[k] = 1
        else
            problem.X_LH[k] = 0
        end
    end
    problem.Z_LH = sum(problem.w[k] * problem.X_LH[k] for k = 1:problem.m)
    return
end


"""
Update Lagrangean multipliers via subgradient method
"""
function update_multipliers!(problem::LagrangeanRelaxedMultiTelescopeTaskingProblem)
    constraint_violation_squared = 0.0
    for k = 1:problem.m
        LHS = sum(sum(problem.A_per_telescope[q][i,k] * problem.Y_LR[q][i] 
                     for i = 1:problem.n_per_telescope[q])
                 for q = 1:problem.s)
        RHS = problem.num_exposure*problem.X_LR[k]
        if LHS < RHS
            constraint_violation_squared += (RHS - LHS)^2
        end
    end

    if constraint_violation_squared < 1e-8
        step = 0.0
    else
        step = problem.Theta * (problem.Z_LR - problem.Z_LH) / constraint_violation_squared
    end

    # update multuplier
    for k = 1:problem.m
        updated_mu = problem.mu[k] + step*(problem.num_exposure*problem.X_LR[k] - 
                                           sum(sum(problem.A_per_telescope[q][i,k] * problem.Y_LR[q][i] 
                                                   for i = 1:problem.n_per_telescope[q])
                                               for q = 1:problem.s))
        problem.mu[k] = max(0.0, updated_mu)
        # problem.mu[k] = updated_mu
    end
    return
end


"""
Solve Multi-telescope tasking problem via Lagrangean Relaxation
"""
function solve_lagrangean_relaxation(
    problem::MultiTelescopeTaskingProblem,
    solver;
    maxiter::Int = 2,
    mu_init::Union{Vector, Nothing} = nothing,
    Theta0::Real = 2.0,
)
    if isnothing(mu_init)
        mu_init = 0.1 * ones(problem.m)
    end

    # initialize Lagrangean relaxation problem
    s = length(problem.n_per_telescope)
    problem_LR = LagrangeanRelaxedMultiTelescopeTaskingProblem(
        mu_init,
        problem.n_per_telescope,
        problem.A_per_telescope,
        problem.T_per_telescope,
        problem.w,
        problem.num_exposure,
        problem.m,
        s,
        # storage for Lagrangean relaxation
        NaN,
        zeros(Int, problem.m),
        Vector{Vector{Int}}(undef, s),
        # storager for Lagrangean heuristic
        NaN,
        zeros(Int, problem.m),
        Vector{Vector{Int}}(undef, s),
        # step size hyperparameter
        Theta0,
    )

    for iter in 1:maxiter
        # solve Lagrangean relaxation problem
        solve_relaxed_problem!(problem_LR, solver)

        # solve Lagrangean heuristic
        heuristic_feasible!(problem_LR)

        # update multipliers
        update_multipliers!(problem_LR)

        # assertions
        # @assert problem_LR.Z_LR >= problem_LR.Z_LH

        gap = (problem_LR.Z_LR - problem_LR.Z_LH)/problem_LR.Z_LR
        @printf("LR iter %d : Z_LR = %1.4f, Z_LH = %1.4f, gap = %1.4f percent \n", iter, problem_LR.Z_LR, problem_LR.Z_LH, 100*gap)
    end

    return problem_LR
end
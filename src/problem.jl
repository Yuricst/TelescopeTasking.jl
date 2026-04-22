"""Struct containing telescope observation scheduling problem information"""


"""
    get_solve_statistics(model::Model)

Get statistics about the model solve into a dictionary
"""
function get_solve_statistics(model::Model)
    stats_dict = Dict(
        "solver_name" => solver_name(model),
        "status" => termination_status(model), #MOI.get(model, MOI.TerminationStatus()),
        "solve_time" => solve_time(model),
        "relative_gap" => relative_gap(model),
        "primal_status" => primal_status(model),
        "dual_status" => dual_status(model),
        "num_variables" => num_variables(model),
    )
    try
        stats_dict["simplex_iterations"] = objective_value(model)
    catch
    end
    try 
        stats_dict["barrier_iterations"] = objective_value(model)
    catch
    end
    return stats_dict
end

struct TelescopeTaskingProblem
    n::Int                      # number of observation arcs
    m::Int                      # number of targets
    w::Vector                   # priority coefficient on each target
    A::Matrix                   # binary n-by-m matrix associating observation arcs to targets
    T::Matrix                   # binary n-by-n matrix for transition feasibility from arc i to j
    num_exposure::Int           # number of exposures required for each target

    function TelescopeTaskingProblem(
        passes::Vector{VisiblePass},
        num_exposure::Int,
        slew_rate::Number;
        w::Union{Vector,Nothing} = nothing,
        buffer_times::Vector = [0.0, 0.0],
    )
        # construct target allocation matrix
        designators = unique([pass.tle.international_designator for pass in passes])
        m = length(designators)
        n = length(passes)
        A = zeros(Int, n, m)
        for (i, pass) in enumerate(passes)
            k = findfirst(x -> x == pass.tle.international_designator, designators)
            A[i,k] = 1
        end
        
        # priority coefficient on each target
        if isnothing(w)
            w = ones(m)
        else
            @assert length(w) == m "Length of w must be equal to number of targets"
        end

        # construct transition feasibility matrix
        T = zeros(Int, n, n)
        for i in 1:n-1
            for j in i:n
                # compute slew time required and slew time allowed
                r_unit_i = sph2cart(vcat(passes[i].azelf_exposure, [1]))
                r_unit_j = sph2cart(vcat(passes[j].azel0_exposure, [1]))
                slew_time = acos(dot(r_unit_i, r_unit_j)) / slew_rate           # in seconds
                slew_allowed = 86400 * (
                    (passes[j].t0_exposure - buffer_times[1]/86400) - (passes[i].tf_exposure + buffer_times[2]/86400)
                )
                if slew_time < slew_allowed
                    T[i,j] = 1
                end
            end
        end
        new(n, m, w, A, T, num_exposure)
    end
end


function Base.show(io::IO, STTP::TelescopeTaskingProblem)
    println(io, "Telescope Scheduling Problem instance")
    println(io, "    Number of passes n  = $(STTP.n)")
    println(io, "    Number of targets m = $(STTP.m)")
end


"""
    solve(problem::TelescopeTaskingProblem, solver; verbose::Bool = true)

Instantiate JuMP model and solve single telescope scheduling problem.
"""
function solve(problem::TelescopeTaskingProblem, solver;
    verbose::Bool = true,
    bias_objective::Bool = false,
    get_model::Bool = false,
)
    if problem.m == 0 
        println("WARNING: no targets detected!")
        return zeros(Int, problem.m), zeros(Int, problem.n), NaN
    end

    # start measuring time
    tstart = time()

    # create model
    model =  Model(solver, add_bridges = false)
    if verbose == false
        set_silent(model)
    end
    @variable(model, Y[1:problem.n], Bin);      # whether observation arc i is selected
    @variable(model, X[1:problem.m], Bin);      # whether target k is observed (sufficiently many times)
    @printf("Created variables; %1.4f sec\n", time() - tstart)

    @constraint(model, sufficient_exposure[k=1:problem.m], 
                sum(problem.A[i,k] * Y[i] for i in 1:problem.n) >= problem.num_exposure * X[k])
    @printf("Created sufficient exposure constraint; %1.4f sec\n", time() - tstart)

    @constraint(model, transition_feasibility[i=1:problem.n-1, j=i+1:problem.n], 
                Y[i] + Y[j] <= 1 + problem.T[i,j])

    @printf("Created transition feasibility constraint; %1.4f sec\n", time() - tstart)
    if bias_objective
        @objective(model, Max, sum(problem.w'X) - 1/problem.n * sum(Y))
    else
        @objective(model, Max, sum(problem.w'X))
    end

    if verbose
        @printf("Model created; elapsed time; %1.4f sec\n", time() - tstart)
    end
    
    if get_model == true
        return model
    end

    # solve problem
    optimize!(model)

    # get BitMatrix X and BitVector Y
    X_val = value.(X) .> 1 - 1e-5;
    Y_val = value.(Y) .> 1 - 1e-5;

    # print info about solution
    if verbose
        termination_status(model)
        @printf("Observed %d out of %d targets, each with %d exposures\n",
            sum(X_val), problem.m, problem.num_exposure)
        @printf("Used %d out of %d passes\n", sum(Y_val), problem.n)
    end
    return X_val, Y_val, get_solve_statistics(model)
end


struct MultiTelescopeTaskingProblem
    n_per_telescope::Vector{Int}        # number of observation arcs per telescope
    m_per_telescope::Vector{Int}        # number of targets per telescope
    n_total::Int                        # total number of observation arcs
    m::Int                              # total number of targets
    A_per_telescope::Vector             # binary n-by-m matrix associating observation arcs to targets per telescope
    T_per_telescope::Vector             # binary n-by-n matrix for transition feasibility from arc i to j per telescopew_per_telescope::Vector             # priority coefficient on each target per telescope
    w::Vector                           # priority coefficient on each target per telescope
    num_exposure::Int                   # number of exposures required for each target
    map_y2Y::Dict                       # index map from individual telescope to stacked Y

    function MultiTelescopeTaskingProblem(
        passes_per_telescope::Vector{Vector{VisiblePass}},
        num_exposure::Int,
        slew_rate::Number;
        w::Union{Vector,Nothing} = nothing,
        buffer_times::Vector = [0.0, 0.0],
    )
        # number of telescopes
        s = length(passes_per_telescope)

        # get weights
        designators = unique([pass.tle.international_designator for pass in vcat(passes_per_telescope...)])
        m = length(designators)

        # priority coefficient on each target
        if isnothing(w)
            w = ones(m)
        else
            @assert length(w) == m "Length of w must be equal to number of targets"
        end
        n_total = sum([length(passes) for passes in passes_per_telescope])
        
        n_per_telescope = Int[]
        m_per_telescope = Int[]
        A_per_telescope = []
        T_per_telescope = []
        for (q,passes) in enumerate(passes_per_telescope)
            # construct target allocation matrix
            n = length(passes)
            A = zeros(Int, n, m)
            observed_names = []
            for (i, pass) in enumerate(passes)
                k = findfirst(x -> x == pass.tle.international_designator, designators)
                A[i,k] = 1
                push!(observed_names, pass.tle.international_designator)
            end
    
            # construct transition feasibility matrix
            T = zeros(Int, n, n)
            for i in 1:n-1
                for j in i:n
                    # compute slew time required and slew time allowed
                    r_unit_i = sph2cart(vcat(passes[i].azelf_exposure, [1]))
                    r_unit_j = sph2cart(vcat(passes[j].azel0_exposure, [1]))
                    slew_time = acos(dot(r_unit_i, r_unit_j)) / slew_rate           # in seconds
                    slew_allowed = 86400 * (
                        (passes[j].t0_exposure - buffer_times[1]/86400) - (passes[i].tf_exposure + buffer_times[2]/86400)
                    )
                    if slew_time < slew_allowed
                        T[i,j] = 1
                    end
                end
            end

            # store into list
            push!(n_per_telescope, n)
            push!(m_per_telescope, length(unique(observed_names)))
            push!(A_per_telescope, A)
            push!(T_per_telescope, T)
        end

        # index map from individual telescope to stacked Y
        i_stacked = 1
        map_y2Y = Dict()
        for q in 1:s
            for i in 1:n_per_telescope[q]
                map_y2Y[(q,i)] = i_stacked
                i_stacked += 1
            end
        end
        new(n_per_telescope, m_per_telescope, n_total, m, A_per_telescope, T_per_telescope, w, num_exposure, map_y2Y)
    end
end


function Base.show(io::IO, MTTP::MultiTelescopeTaskingProblem)
    println(io, "Telescopes Scheduling Problem instance")
    println(io, "    Number of telescopes s = $(length(MTTP.n_per_telescope))")
    for q in 1:length(MTTP.n_per_telescope)
        println(io, "    Telescope $q")
        println(io, "        Number of passes n  = $(MTTP.n_per_telescope[q])")
        println(io, "        Number of targets m = $(MTTP.m_per_telescope[q])")
        println(io, "        size(A) = $(size(MTTP.A_per_telescope[q]))")
        println(io, "        size(T) = $(size(MTTP.T_per_telescope[q]))")
    end
end


"""
    decompose_multitelescope_Y(n_per_telescope::Vector{Int}, map_y2Y::Dict, Y::Union{Vector,BitVector})

Decompose optimal solution into per-telescope solutions
"""
function decompose_multitelescope_Y(n_per_telescope::Vector{Int}, map_y2Y::Dict, Y::Union{Vector,BitVector})
    Y_per_telescope = Vector[]
    for q in 1:length(n_per_telescope)
        Y_per_telescope_q = zeros(Int, n_per_telescope[q])
        for i in 1:n_per_telescope[q]
            Y_per_telescope_q[i] = Y[map_y2Y[(q,i)]]
        end
        push!(Y_per_telescope, Y_per_telescope_q)
    end
    return Y_per_telescope
end


"""
    decompose_multitelescope_Y(problem::MultiTelescopeTaskingProblem, Y::Union{Vector,BitVector})

Decompose optimal solution into per-telescope solutions
"""
function decompose_multitelescope_Y(problem::MultiTelescopeTaskingProblem, Y::Union{Vector,BitVector})
    return decompose_multitelescope_Y(problem.n_per_telescope, problem.map_y2Y, Y)
end


"""
    solve(problem::MultiTelescopeTaskingProblem, solver; verbose::Bool = true)

Instantiate JuMP model and solve multi telescope scheduling problem.
"""
function solve(
    problem::MultiTelescopeTaskingProblem,
    solver;
    verbose::Bool = true,
    bias_objective::Bool = false,
    get_model::Bool = false,
)
    # start measuring time
    tstart = time()

    # create model
    model = Model(solver; add_bridges = false)
    if verbose == false
        set_silent(model)
    end
    #N = sum(problem.n_per_telescope)
    @variable(model, Y[1:problem.n_total], Bin);      # whether observation arc i is selected
    @variable(model, X[1:problem.m], Bin);            # whether target k is observed (sufficiently many times)
    @printf("Created variables; %1.4f sec\n", time() - tstart)

    # sufficient exposure constraint
    s = length(problem.n_per_telescope)
    @constraint(model,
                sufficient_exposure[k=1:problem.m], 
                sum(
                    sum(
                        problem.A_per_telescope[q][i,k] * Y[problem.map_y2Y[(q,i)]]
                        for i in 1:problem.n_per_telescope[q]
                    )
                for q in 1:s) >= problem.num_exposure * X[k]
    )
    @printf("Created target allocation constraints; %1.4f sec\n", time() - tstart)
    
    # transition feasibility constraints
    @showprogress for (q,T) in enumerate(problem.T_per_telescope)
        @constraint(model,
                   [i = 1:problem.n_per_telescope[q]-1, j = i+1:problem.n_per_telescope[q]], 
                    Y[problem.map_y2Y[(q,i)]] + Y[problem.map_y2Y[(q,j)]] <= 1 + T[i,j])
    end
    @printf("Created transition feasibility constraints; %1.4f sec\n", time() - tstart)

    if bias_objective == true
        @objective(model, Max, sum(problem.w'X) - 1/problem.n_total * sum(Y))
    else
        @objective(model, Max, sum(problem.w'X))
    end

    if verbose
        @printf("Model created; elapsed time; %1.4f sec\n", time() - tstart)
    end

    if get_model == true
        return model
    end

    # solve problem
    optimize!(model)

    # get BitMatrix X and BitVector Y
    X_val = value.(X) .> 1 - 1e-5;
    Y_val = value.(Y) .> 1 - 1e-5;

    # decompose solution
    Y_per_telescope = decompose_multitelescope_Y(problem, Y_val)

    # print info about solution
    if verbose
        termination_status(model)
        @printf("Observed %d out of %d targets, each with %d exposures\n",
            sum(X_val), problem.m, problem.num_exposure)
        @printf("Used %d out of %d passes\n", sum(Y_val), problem.n_total)
        for q = 1:s
            @printf("Telescope %d: %d out of %d passes used\n", 
                q, sum(Y_per_telescope[q]), problem.n_per_telescope[q])
        end
    end
    return X_val, Y_val, Y_per_telescope, get_solve_statistics(model)
end


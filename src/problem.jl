"""Struct containing telescope observation scheduling problem information"""


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
            j = findfirst(x -> x == pass.tle.international_designator, designators)
            A[i,j] = 1
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


function Base.show(io::IO, tsp::TelescopeTaskingProblem)
    println(io, "Telescope Scheduling Problem instance")
    println(io, "    Number of observation arcs n = $(tsp.n)")
    println(io, "    Number of targets m = $(tsp.m)")
end


"""
    solve!(problem::TelescopeTaskingProblem, solver; verbose::Bool = true)

Instantiate JuMP model and solve telescope scheduling problem.
"""
function solve!(problem::TelescopeTaskingProblem, solver; verbose::Bool = true)
    if problem.m == 0 
        println("WARNING: no targets detected!")
        return zeros(Int, problem.m), zeros(Int, problem.n), NaN
    end

    # start measuring time
    tstart = time()

    # create model
    model = Model(solver)
    @variable(model, Y[1:problem.n], Bin);      # whether observation arc i is selected
    @variable(model, X[1:problem.m], Bin);      # whether target k is observed (sufficiently many times)
    @printf("Created variables; %1.4f sec\n", time() - tstart)

    @constraint(model, sufficient_exposure[k=1:problem.m], 
                sum(problem.A[i,k] * Y[i] for i in 1:problem.n) >= problem.num_exposure * X[k])
    @printf("Created sufficient exposure constraint; %1.4f sec\n", time() - tstart)

    @constraint(model, transition_feasibility[i=1:problem.n-1, j=i+1:problem.n], 
                Y[i] + Y[j] <= 1 + problem.T[i,j])

    @printf("Created transition feasibility constraint; %1.4f sec\n", time() - tstart)
    @objective(model, Max, sum(problem.w'X) - 1/problem.n * sum(Y))

    if verbose
        @printf("Model created; elapsed time; %1.4f sec\n", time() - tstart)
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
    return X_val, Y_val, MOI.get(model, MOI.TerminationStatus())
end
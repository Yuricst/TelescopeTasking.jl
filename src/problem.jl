"""Struct containing telescope observation scheduling problem information"""


struct TelescopeSchedulingProblem
    n::Int                      # number of observation arcs
    m::Int                      # number of targets
    A::Matrix                   # binary n-by-m matrix associating observation arcs to targets
    T::Matrix                   # binary n-by-n matrix for transition feasibility from arc i to j
    num_exposure::Int           # number of exposures required for each target

    function TelescopeSchedulingProblem(
        passes::Vector{VisiblePass},
        num_exposure::Int,
        slew_rate::Number;
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
        new(n, m, A, T, num_exposure)
    end
end


function Base.show(io::IO, tsp::TelescopeSchedulingProblem)
    println(io, "Telescope Scheduling Problem instance")
    println(io, "    Number of observation arcs n = $(tsp.n)")
    println(io, "    Number of targets m = $(tsp.m)")
end


"""
    solve!(problem::TelescopeSchedulingProblem, solver)

Instantiate JuMP model and solve telescope scheduling problem.
"""
function solve!(problem::TelescopeSchedulingProblem, solver; verbose::Bool = true)

    # create model
    model = Model(solver)
    @variable(model, Y[1:problem.n], Bin);      # whether observation arc i is selected
    @variable(model, X[1:problem.m], Bin);      # whether target k is observed (sufficiently many times)
    
    @constraint(model, sufficient_exposure[k=1:problem.m], 
                sum(problem.A[i,k] * Y[i] for i in 1:problem.n) >= problem.num_exposure * X[k])

    @constraint(model, transition_feasibility[i=1:problem.n-1, j=i+1:problem.n], 
                Y[i] + Y[j] <= 1 + problem.T[i,j])

    @objective(model, Max, sum(X) - 1/problem.m * sum(Y))
    optimize!(model)

    # get BitMatrix X and BitVector Y
    X_val = value.(X) .> 1 - 1e-5;
    Y_val = value.(Y) .> 1 - 1e-5;

    # print info about solution
    if verbose
        termination_status(model)
        @printf("Number of observed targets: %d out of %d in FOV\n", sum(value.(X)), problem.m)
    end
    return X_val, Y_val
end
"""Run all tests"""

@time begin
    @time "Passes IO" include("test_passes_io.jl")
    @time "Solve single telescope problem" include("test_solve_single.jl")
    @time "Solve multi telescope problem" include("test_solve_MTTP.jl")
end

println("Done!")
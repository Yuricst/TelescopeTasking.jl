"""Run all tests"""

@time begin
    @time "Passes IO" include("test_passes_io.jl")
end

println("Done!")
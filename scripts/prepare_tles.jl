"""Prepare TLEs to be used in specific configuration"""

using Printf: @printf
using SatelliteToolboxTle

include(joinpath(@__DIR__, "../src/TelescopeTasking.jl"))

function tles_config_A()
    println("\n ********* Generating TLEs for config A ********* ")
    
    path_to_tles = joinpath(@__DIR__, "..", "data", "tles", "active.txt")
    tles_str = read(path_to_tles, String)

    # convert to TLE objects
    tles = read_tles(tles_str)
    @printf("%d objects from active.txt loaded\n", length(tles))

    # filter them
    names_include = ["GLOBALSTAR", "IRIDIUM"]
    tles = TelescopeTasking.filter(tles, names_include = names_include)
    @printf("%d objects after filtering for ", length(tles))
    println(names_include)

    # save to file tles used
    TelescopeTasking.tles_to_file(tles, joinpath(@__DIR__, "../data/tles/AAS25targetA.txt"))
end


function tles_config_B()
    println("\n ********* Generating TLEs for config B ********* ")
    
    tles = SatelliteToolboxTle.TLE[]
    files = ["cosmos-1408-debris.txt", "cosmos-2251-debris.txt", "fengyun-1c-debris.txt", "iridium-33-debris.txt"]
    for file in files
        _path_to_tles = joinpath(@__DIR__, "..", "data", "tles", file)
        _tles_str = read(_path_to_tles, String)
        append!(tles, read_tles(_tles_str))
    end
    @printf("%d objects loaded from files ", length(tles))
    println(files)

    # save to file tles used
    TelescopeTasking.tles_to_file(tles, joinpath(@__DIR__, "../data/tles/AAS25targetB.txt"))
end


function tles_config_C()
    println("\n ********* Generating TLEs for config C ********* ")
    
    path_to_tles = joinpath(@__DIR__, "..", "data", "tles", "active.txt")
    tles_str = read(path_to_tles, String)

    # convert to TLE objects
    tles = read_tles(tles_str)
    @printf("%d objects from active.txt loaded\n", length(tles))

    # filter them
    names_include = ["STARLINK",]
    tles = TelescopeTasking.filter(tles, names_include = names_include)
    @printf("%d objects after filtering for ", length(tles))
    println(names_include)

    # save to file tles used
    TelescopeTasking.tles_to_file(tles, joinpath(@__DIR__, "../data/tles/AAS25targetC.txt"))
    println("Done!")
end

tles_config_A()
tles_config_B()
tles_config_C()
println("Done!")
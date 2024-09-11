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

    # print latest epoch
    jd0_ref = maximum([tle_epoch(tle) for tle in tles])
    @printf("Latest epoch JD: %1.6f\n", jd0_ref)
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

    # print latest epoch
    jd0_ref = maximum([tle_epoch(tle) for tle in tles])
    @printf("Latest epoch JD: %1.6f\n", jd0_ref)
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
    
    # print latest epoch
    jd0_ref = maximum([tle_epoch(tle) for tle in tles])
    @printf("Latest epoch JD: %1.6f\n", jd0_ref)
end


function tles_config_D()
    println("\n ********* Generating TLEs for config D ********* ")
    
    tles = SatelliteToolboxTle.TLE[]
    files = ["cosmos-2251-debris.txt", "iridium-33-debris.txt"]
    for file in files
        _path_to_tles = joinpath(@__DIR__, "..", "data", "tles", file)
        _tles_str = read(_path_to_tles, String)
        append!(tles, read_tles(_tles_str))
    end
    @printf("%d objects loaded from files ", length(tles))
    println(files)

    # save to file tles used
    TelescopeTasking.tles_to_file(tles, joinpath(@__DIR__, "../data/tles/AAS25targetD.txt"))

    # print latest epoch
    jd0_ref = maximum([tle_epoch(tle) for tle in tles])
    @printf("Latest epoch JD: %1.6f\n", jd0_ref)
end


function tles_config_S(;subname::String)
    println("\n ********* Generating TLEs for config S$subname ********* ")
    
    path_to_tles = joinpath(@__DIR__, "..", "data", "tles", "active.txt")
    tles_str = read(path_to_tles, String)

    # convert to TLE objects
    tles = read_tles(tles_str)
    @printf("%d objects from active.txt loaded\n", length(tles))

    # filter them
    names_include = ["STARLINK",]
    tles = TelescopeTasking.filter(tles, names_include = names_include)
    if subname == "1"
        tles = filter(tle -> tle.inclination < 50, tles)
    elseif subname == "2"
        tles = filter(tle -> 50 < tle.inclination < 60, tles)
    elseif subname == "3"
        tles = filter(tle -> 60 < tle.inclination < 80, tles)
    elseif subname == "4"
        tles = filter(tle -> 80 < tle.inclination, tles)
    end
    @printf("%d objects after filtering for ", length(tles))
    println(names_include)

    # save to file tles used
    TelescopeTasking.tles_to_file(tles, joinpath(@__DIR__, "../data/tles/AAS25targetS$subname.txt"))
    
    # print latest epoch
    jd0_ref = maximum([tle_epoch(tle) for tle in tles])
    @printf("Latest epoch JD: %1.6f\n", jd0_ref)
end


# tles_config_A()
# tles_config_B()
# tles_config_C()
# tles_config_D()
tles_config_S(subname = "0")
tles_config_S(subname = "1")
tles_config_S(subname = "2")
tles_config_S(subname = "3")
tles_config_S(subname = "4")
println("Done!")
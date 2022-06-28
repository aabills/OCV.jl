using Test
@testset "test ocv speed and allocations" begin include("ocv_speedtests.jl") end
@testset "test inversion to find filling fractions" begin include("ocv_inversiontest.jl") end

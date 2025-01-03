using LatticeGeometry
using Test

@testset "LatticeGeometry.jl" begin
    # Write your tests here.
    include("test_UnitCells.jl")
    include("test_CrystalGeometer.jl")
    include("test_CellIndexing.jl")
    include("test_LatticeIndexing.jl")
end

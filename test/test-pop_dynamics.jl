@testset "Population dynamics" begin
    @test sum(P1) ≈ K
    @test all(vecstate(P1) .== vecstate(eqdist_ap))
    @test all(vecstate(Phalf1) .≥ vecstate(Phalf))
    @test sum(Phalf1) ≈ 51.25

    @test P1 isa PopState
    @test P1_pt isa PopState
    @test P1_uspt isa PopState
    @test P1_mspt isa PopState
end


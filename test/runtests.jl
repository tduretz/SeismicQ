using SeismicQ, Test

@testset "First series of tests" verbose=true begin
    @test 1==1
    @test Ricker(0.0, 0.0, 0.0) ≈ 1.0
end
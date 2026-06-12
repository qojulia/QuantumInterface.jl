using QuantumInterface
using Test
using JET

using JET: MethodErrorReport

@testset "JET checks" begin
    rep = report_package(QuantumInterface; target_modules=(QuantumInterface,))
    @show rep
    reports = JET.get_reports(rep)
    non_method_reports = filter(report -> !(report isa MethodErrorReport), reports)
    @test isempty(non_method_reports)
    @test_broken length(JET.get_reports(rep)) == 0
end

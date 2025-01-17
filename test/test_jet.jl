using QuantumInterface
using Test
using JET

using JET: ReportPass, DefinitionAnalysisPass, BasicPass, InferenceErrorReport, UncaughtExceptionReport, MethodErrorReport

# We define abstract methods, so it is unsurprising that there are no matching methods unless
# we import other libraries.
# TODO find a more fine-grained way to check these.
struct NoMatchingMethodIsOK <: ReportPass end

# ignores `MethodErrorReport` analyzed by `JETAnalyzer`
(::NoMatchingMethodIsOK)(::Type{MethodErrorReport}, @nospecialize(_...)) = return

function (::NoMatchingMethodIsOK)(report_type::Type{<:InferenceErrorReport}, @nospecialize(args...))
    DefinitionAnalysisPass()(report_type, args...) # Using it instead of BasePass, see https://github.com/aviatesk/JET.jl/pull/532
end

@testset "JET checks" begin
    rep = report_package("QuantumInterface";
        report_pass=NoMatchingMethodIsOK(),
        ignored_modules=(
            #AnyFrameModule(...),
        )
    )
    @show rep
    @test length(JET.get_reports(rep)) <= 3
    @test_broken length(JET.get_reports(rep)) == 0
end

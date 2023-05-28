using QuantumInterface
using JET

using JET: ReportPass, BasicPass, InferenceErrorReport, UncaughtExceptionReport, MethodErrorReport

# Custom report pass that ignores `UncaughtExceptionReport`
# Too coarse currently, but it serves to ignore the various
# "may throw" messages for runtime errors we raise on purpose
# (mostly on malformed user input)
struct MayThrowIsOk <: ReportPass end

# We define abstract methods, so it is unsurprising that there are no matching methods unless
# we import other libraries.
# TODO find a more fine-grained way to check these.
struct NoMatchingMethodIsOK <: ReportPass end

# ignores `UncaughtExceptionReport` analyzed by `JETAnalyzer`
(::MayThrowIsOk)(::Type{UncaughtExceptionReport}, @nospecialize(_...)) = return

# forward to `BasicPass` for everything else
function (::MayThrowIsOk)(report_type::Type{<:InferenceErrorReport}, @nospecialize(args...))
    BasicPass()(report_type, args...)
end

# ignores `UncaughtExceptionReport` and `MethodErrorReport` analyzed by `JETAnalyzer`
(::NoMatchingMethodIsOK)(::Type{MethodErrorReport}, @nospecialize(_...)) = return

# forward to `MayThrowIsOk` for everything else
function (::NoMatchingMethodIsOK)(report_type::Type{<:InferenceErrorReport}, @nospecialize(args...))
    MayThrowIsOk()(report_type, args...)
end

@testset "JET checks" begin
    rep = report_package("QuantumInterface";
        report_pass=MayThrowIsOk(),
        ignored_modules=(
            #AnyFrameModule(...),
        )
    )
    @show rep
    @test length(JET.get_reports(rep)) <= 4

    rep = report_package("QuantumInterface";
        report_pass=NoMatchingMethodIsOK(),
        ignored_modules=(
            #AnyFrameModule(...),
        )
    )
    @show rep
    @test length(JET.get_reports(rep)) == 0
end

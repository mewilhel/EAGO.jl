using BenchmarkTools, MINLPLib, BenchmarkProfiles

T = 10 * rand(25, 3)
performance_profile(T, ["SCIP", "EAGO", "APLINE"])

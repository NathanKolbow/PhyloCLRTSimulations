
"""
Likelihood components of the quartet (SNaQ) model, in the conventions the PhyloCLRT tests
expect (the same as the new models' `logf`, `sensitivity` and `variability`):
- `logf` and `grad` are SUMMED over the `length(gts)` loci. `computeSNaQscore!` and
  `computegradient` work on CF proportions, i.e. they return the per-locus average.
- `sens` is the Hessian itself (negative definite at a maximum); the tests use H = -sens.

`computeHessian` and `compute_Jacobian` share a 1/#quartets factor, which cancels in every
test statistic.
"""
function quartet_likelihood_components(net::HybridNetwork, gts::Vector{HybridNetwork}, ocfs::Matrix{Float64}, ρ::Float64)
	eqns = SNaQ.findquartetequations(net)[1];
	var = PhyloCLRT.compute_Jacobian(net, eqns, gts, ρ)
	sens = PhyloCLRT.computeHessian(net, eqns, ocfs, ρ)
	grad = length(gts) .* SNaQ.computegradient(net, ocfs, ρ)
	logf = length(gts) * SNaQ.computeSNaQscore!(net, ocfs, ρ)
	return logf, grad, sens, var
end

function old_model_comptypes(net0::HybridNetwork, net1::HybridNetwork, gts::Vector{HybridNetwork}, ρ::Float64)::CompTypes
	ocfs = gts2CFs(gts)
	logf1, grad1, sens1, var1 = quartet_likelihood_components(net1, gts, ocfs, ρ)
	logf0, _, sens0, var0 = quartet_likelihood_components(net0, gts, ocfs, ρ)
	return (logf0, sens0, var0, logf1, grad1, sens1, var1)
end

run_test_old_model(net0::HybridNetwork, net1::HybridNetwork, gts::Vector{HybridNetwork}, test::String, ρ::Float64=1.0) =
	run_test(net1, test, old_model_comptypes(net0, net1, gts, ρ))

function run_tests_old_model(net0::HybridNetwork, net1::HybridNetwork, gts::Vector{HybridNetwork}, tests::Vector{String}, ρ::Float64=1.0)::Vector{Float64}
	comps = old_model_comptypes(net0, net1, gts, ρ)
	return Float64[run_test(net1, test, comps) for test in tests]
end

function quartetCLICstatistic(net::HybridNetwork, gts::Vector{HybridNetwork}, ρ::Float64=1.0)
	logf, _, sens, var = quartet_likelihood_components(net, gts, gts2CFs(gts), ρ)
	return CLICstatistic(var, sens, logf)
end

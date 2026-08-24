
function run_test_old_model(net0::HybridNetwork, net1::HybridNetwork, gts::Vector{HybridNetwork}, test::String, ρ::Float64=1.0)
	ocfs = gts2CFs(gts)

	eqns1 = SNaQ.findquartetequations(net1)[1];
	var1 = PhyloCLRT.compute_Jacobian(net1, eqns1, gts, ρ)
	sens1 = .-PhyloCLRT.computeHessian(net1, eqns1, ocfs, ρ)
	grad1 = SNaQ.computegradient(net1, ocfs, ρ)
	logf1 = SNaQ.computeSNaQscore!(net1, ocfs, ρ)
	
	eqns0 = SNaQ.findquartetequations(net0)[1];
	var0 = PhyloCLRT.compute_Jacobian(net0, eqns0, gts, ρ)
	sens0 = .-PhyloCLRT.computeHessian(net0, eqns0, ocfs, ρ)
	logf0 = SNaQ.computeSNaQscore!(net0, ocfs, ρ)

	return run_test(net1, test, CompTypes([
		logf0, sens0, var0, logf1, grad1, sens1, var1
	]))
end

function run_tests_old_model(net0::HybridNetwork, net1::HybridNetwork, gts::Vector{HybridNetwork}, tests::Vector{String}, ρ::Float64=1.0)::Vector{Float64}
	ocfs = gts2CFs(gts)

	eqns1 = SNaQ.findquartetequations(net1)[1];
	var1 = PhyloCLRT.compute_Jacobian(net1, eqns1, gts, ρ)
	sens1 = .-PhyloCLRT.computeHessian(net1, eqns1, ocfs, ρ)
	grad1 = SNaQ.computegradient(net1, ocfs, ρ)
	logf1 = computeSNaQscore!(net1, ocfs, ρ)
	
	eqns0 = SNaQ.findquartetequations(net0)[1];
	var0 = PhyloCLRT.compute_Jacobian(net0, eqns0, gts, ρ)
	sens0 = .-PhyloCLRT.computeHessian(net0, eqns0, ocfs, ρ)
	logf0 = computeSNaQscore!(net0, ocfs, ρ)

	return Float64[run_test(net1, test, CompTypes([logf0, sens0, var0, logf1, grad1, sens1, var1])) for test in tests]
end

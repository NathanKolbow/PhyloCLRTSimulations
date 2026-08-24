const CompTypes = Tuple{Float64, Matrix{Float64}, Matrix{Float64}, Float64, Vector{Float64}, Matrix{Float64}, Matrix{Float64}};
SMOOTH_OPTIONS = ["smooth", "rough", "smooth_sensitivity", "smooth_variability", "smooth_sensvar"];

function gather_likelihood_components(net::HybridNetwork, gts::Vector{HybridNetwork}, model::String, eps::Float64, smoothness::String=ϵ == 0.0 ? "rough" : "smooth")
	eps >= 0 || error("eps must be at least 0 (eps=$eps)")
	(eps == 0.0 && smoothness == "rough") || eps > 0 || error("smoothness must be \"rough\" for eps=0.0 (smoothness=$smoothness)")

	ds = if model == "joint"
		gather_mixture_distributions(net)
	elseif model == "marginal"
		MyDistMarginalized.(gather_mixture_distributions(net))
	else
		error("model $model not recognized")
	end

	ts = if model == "joint"
		t_distribution(gts)
	elseif model == "marginal"
		marginalize_tD(t_distribution(gts))
	else
		error("model $model not recognized")
	end

	return if eps == 0.0
		logf(ds, ts),
		∇logf(ds, ts),
		sensitivity(ds, gts),
		variability(ds, gts)
	elseif smoothness == "smooth"
		smooth_logf(ds, ts, eps),
		smooth_∇logf(ds, ts, eps),
		smooth_sensitivity(ds, gts, eps),
		smooth_variability(ds, gts, eps)
	elseif smoothness == "smooth_sensitivity"
		logf(ds, ts),
		∇logf(ds, ts),
		smooth_sensitivity(ds, gts, eps),
		variability(ds, gts)
	elseif smoothness == "smooth_variability"
		logf(ds, ts),
		∇logf(ds, ts),
		sensitivity(ds, gts),
		smooth_variability(ds, gts, eps)
	elseif smoothness == "smooth_sensvar"
		logf(ds, ts),
		∇logf(ds, ts),
		smooth_sensitivity(ds, gts, eps),
		smooth_variability(ds, gts, eps)
	else
		error("smoothness \"$smoothness\" not recognized")
	end
end

function gather_likelihood_components(
	distH0::Union{Vector{MyDistConcat}, Vector{MyDistMarginalized}, MyDistMarginalizedAllQuartets},
	distH1::Union{Vector{MyDistConcat}, Vector{MyDistMarginalized}, MyDistMarginalizedAllQuartets},
	gts::Vector{HybridNetwork},
	smoothness::String,
	eps::Float64)::CompTypes

	smoothness in SMOOTH_OPTIONS || error("smoothness option \"$(smoothness)\" not recognized.")
	typeof(distH0) == typeof(distH1) || error("distH0 is $(typeof(distH0)) but distH1 is $(typeof(distH1))")

	ts = t_distribution(gts);
	if typeof(distH0) <: Vector{MyDistMarginalized}
		ts = marginalize_tD(ts)
	elseif typeof(distH0) <: MyDistMarginalizedAllQuartets
		ts = reduce(vcat, ts)
	end

	return if smoothness == "smooth"
		smooth_logf(distH0, ts, eps),
		smooth_sensitivity(distH0, gts, eps),
		smooth_variability(distH0, gts, eps),
		smooth_logf(distH1, ts, eps),
		smooth_∇logf(distH1, ts, eps),
		smooth_sensitivity(distH1, gts, eps),
		smooth_variability(distH1, gts, eps)
	elseif smoothness == "rough"
		logf(distH0, ts),
		sensitivity(distH0, gts),
		variability(distH0, gts),
		logf(distH1, ts),
		∇logf(distH1, ts),
		sensitivity(distH1, gts),
		variability(distH1, gts)
	elseif smoothness == "smooth_sensitivity"
		logf(distH0, ts),
		smooth_sensitivity(distH0, gts, eps),
		variability(distH0, gts),
		logf(distH1, ts),
		∇logf(distH1, ts),
		smooth_sensitivity(distH1, gts, eps),
		variability(distH1, gts)
	elseif smoothness == "smooth_variability"
		logf(distH0, ts),
		sensitivity(distH0, gts),
		smooth_variability(distH0, gts, eps),
		logf(distH1, ts),
		∇logf(distH1, ts),
		sensitivity(distH1, gts),
		smooth_variability(distH1, gts, eps)
	elseif smoothness == "smooth_sensvar"
		logf(distH0, ts),
		smooth_sensitivity(distH0, gts, eps),
		smooth_variability(distH0, gts, eps),
		logf(distH1, ts),
		∇logf(distH1, ts),
		smooth_sensitivity(distH1, gts, eps),
		smooth_variability(distH1, gts, eps)
	else
		error("\"$(smoothness)\" not a recognized smoothness value.")
	end
end

function select_likelihood_components(
	smoothness::String,
	rlogf0::Float64, rsens0::Matrix{Float64}, rvar0::Matrix{Float64},
	rlogf1::Float64, rgrad1::Vector{Float64}, rsens1::Matrix{Float64}, rvar1::Matrix{Float64},
	slogf0::Float64, ssens0::Matrix{Float64}, svar0::Matrix{Float64},
	slogf1::Float64, sgrad1::Vector{Float64}, ssens1::Matrix{Float64}, svar1::Matrix{Float64}
)::CompTypes
	return if smoothness == "smooth"
		slogf0, ssens0, svar0, slogf1, sgrad1, ssens1, svar1
	elseif smoothness == "rough"
		rlogf0, rsens0, rvar0, rlogf1, rgrad1, rsens1, rvar1
	elseif smoothness == "smooth_sensitivity"
		rlogf0, ssens0, rvar0, rlogf1, rgrad1, ssens1, rvar1
	elseif smoothness == "smooth_variability"
		rlogf0, rsens0, svar0, rlogf1, rgrad1, rsens1, svar1
	elseif smoothness == "smooth_sensvar"
		rlogf0, ssens0, svar0, rlogf1, rgrad1, ssens1, svar1
	else
		error("$(smoothness) not a recognized smoothness value.")
	end
end

function get_snaq_param_errors(N::HybridNetwork, dcf::DataCF)::Float64
	N = SNaQ.deepcopynetwork(N);
	SNaQ.semidirectnetwork!(N)
	trueparams = SNaQ.gatherparams(N);
	fitnumericalparameters!(N, dcf, 1.0)
	return sum(abs.(SNaQ.gatherparams(N) .- trueparams))
end

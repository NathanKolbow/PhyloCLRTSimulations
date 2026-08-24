include("../helpers/matrix-selection-helpers.jl")


function get_dists(tre::HybridNetwork, net::HybridNetwork, model::String)
	return if model == "joint"
		gather_mixture_distributions(tre), gather_mixture_distributions(net)
	elseif model == "marginal"
		MyDistMarginalized.(gather_mixture_distributions(tre)), MyDistMarginalized.(gather_mixture_distributions(net))
	elseif model == "supermarginal"
		MyDistMarginalizedAllQuartets(MyDistMarginalized.(gather_mixture_distributions(tre))), MyDistMarginalizedAllQuartets(MyDistMarginalized.(gather_mixture_distributions(net)))
	else
		error("Model type $(model) not recognized.")
	end
end

"""
Returns a p-value for proper tests, or the difference in IC values for ICs.
"""
function effect_size_test(T::HybridNetwork, N::HybridNetwork, test::String, lik_components::CompTypes)::Float64
	T.numhybrids == 0 || error("T is FIRST argument.")
	N.numhybrids > 0 || error("N is SECOND argument.")

	logf0, sens0, var0, logf1, grad1, sens1, var1 = lik_components;
	return if test == "CLIC"
		CLICstatistic(var1, sens1, logf1) - CLICstatistic(var0, sens0, logf0)
	elseif test == "cw"
		cw(var1, sens1, logf1, logf0)[2]
	elseif test == "cwP"
		cwP(N, var1, sens1, logf1, logf0)[2]
	elseif test == "cLRI"
		cLRI(var1, sens1, grad1, logf1, logf0)[2]
	elseif test == "cLR1"
		cLR1(var1, sens1, logf1, logf0)
	elseif test == "cLR2"
		cLR2(var1, sens1, logf1, logf0)
	elseif test == "cLR"
		cLR(var1, sens1, logf1, logf0)
	else
		error("test $(test) not recognized")
	end
end

function has(d::DataFrame, ϵ::Float64, model::String, γ::Float64, ngt::Int64, alltests::Vector{String})::Bool
	r = filter(r -> r.eps == ϵ && r.model == model && r.gamma == γ && r.ngt == ngt, d)
	return all(t -> t in r.test, alltests)
end

function has(d::DataFrame, ϵ::Float64, model::String, γ::Float64, test::String)::Bool
	r = filter(r -> r.eps == ϵ && r.model == model && r.gamma == γ, d)
	return test in r.test
end

function hasall(d::DataFrame, ϵ::Float64, model::String, ngt::Int64, alltests::Vector{String}, nreps::Int64)::Bool
	r = filter(r -> r.eps == ϵ && r.model == model && r.ngt == ngt, d)
	return all(t -> sum(r.test .== t) >= nreps, alltests)
end

function setγs(net::HybridNetwork, γ::Float64)
	for H in net.hybrid
		getparentedge(H).gamma = 1.0 - γ
		getparentedgeminor(H).gamma = γ
	end
end

function run_test(H1net::HybridNetwork, test::String, lik_components::CompTypes)::Float64
	logf0, sens0, var0, logf1, grad1, sens1, var1 = lik_components;
	return if test == "CLIC"
		CLICstatistic(var1, sens1, logf1) - CLICstatistic(var0, sens0, logf0)
	elseif test == "cw"
		cw(var1, sens1, logf1, logf0)[2]
	elseif test == "cwP"
		cwP(H1net, var1, sens1, logf1, logf0)[2]
	elseif test == "cLRI"
		cLRI(var1, sens1, grad1, logf1, logf0)[2]
	elseif test == "cLR1"
		cLR1(var1, sens1, logf1, logf0)
	elseif test == "cLR2"
		cLR2(var1, sens1, logf1, logf0)
	elseif test == "cLR"
		cLR(var1, sens1, logf1, logf0)
	else
		error("test $(test) not recognized")
	end
end

function optimize_given_model(net::HybridNetwork, gts::Vector{HybridNetwork}, model::String, ϵ::Float64; verbose::Bool=false, rootmaxeval=100, finalmaxeval=1000)::HybridNetwork
	model in ["joint", "marginal"] || error("Only joint and marginal models accepted (model=$model).")
	net = SNaQ.deepcopynetwork(net);

	if model == "joint"
		net = optimize_root_placement(net, gts, ϵ; maxeval=rootmaxeval, verbose=verbose)
		LL = optimize_parameters∇!(net, t_distribution(gts), ϵ; maxeval=finalmaxeval, verbose=verbose)
		SNaQscore!(net, LL)
	elseif model == "marginal"
		net = mdm_optimize_root_placement(net, gts, ϵ; maxeval=rootmaxeval, verbose=verbose)
		LL = mdm_optimize_parameters∇!(net, t_distribution(gts), ϵ; maxeval=finalmaxeval, verbose=verbose)
		SNaQscore!(net, LL)
	else
		error("$model not recognized model.")
	end
	return net
end

function snaqnetsearch(
	truenet::HybridNetwork, dcf::DataCF, hmax::Int64, ρ::Float64=1.0;
	Nfail::Int64=10, runs::Int64=100
)::HybridNetwork
	net = SNaQ.deepcopynetwork(truenet);
	if hmax == net.numhybrids
		fitnumericalparameters!(net, dcf, ρ)
	elseif hmax > net.numhybrids
		net = snaq!(net, dcf; hmax=hmax, runs=runs, Nfail=Nfail, filename="", ρ=ρ)
	elseif hmax == 0
		net = majortree(net)
		fitnumericalparameters!(net, dcf, ρ)
	elseif hmax < net.numhybrids
		for j = 1:(net.numhybrids - hmax)
			bestproposal = nothing
			for j in eachindex(net.hybrid)
				prop = SNaQ.deepcopynetwork(net)
				SNaQ.removehybrid!(prop, prop.hybrid[j], true)
				SNaQ.semidirectnetwork!(prop)
				SNaQ.fitnumericalparameters!(prop, dcf, ρ)
				if isnothing(bestproposal) || SNaQscore(prop) > SNaQscore(bestproposal)
					bestproposal = prop
				end
			end
			net = bestproposal
		end
		SNaQ.fitnumericalparameters!(net, dcf, ρ)
	end
	return net
end
using SNaQ
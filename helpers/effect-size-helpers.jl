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
	return if test == "cw"
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

"""
Returns 2-tuple with (abs sum of internal edge length errors, abs gamma error) between two
networks with the same topology and 1 hybrid.

Edges are paired with root-free keys (see `rootfreeedgekeys`), so the result does not depend on
where either network is rooted, on whether it is semidirected, or on which hybrid parent edge
is labeled major. Pendant edges are not estimated by either model and are ignored.
"""
function absparamerrors(truenet::HybridNetwork, estnet::HybridNetwork)::NTuple{2, Float64}
	truenet.numhybrids == 1 && estnet.numhybrids == 1 || error("Function written for only 1 hybrid")
	truechains = rootfreeedgekeys(truenet)
	estchains = rootfreeedgekeys(estnet)
	Set(keys(truechains)) == Set(keys(estchains)) || error("truenet and estnet have different topologies")

	absterror = 0.0
	for (k, c) in truechains
		c.pendant && continue
		absterror += abs(c.length - estchains[k].length)
	end

	# γ of the true minor hybrid edge vs. γ of the same edge in `estnet`
	trueminor = getparentedgeminor(truenet.hybrid[1])
	k = only(k for (k, c) in truechains if any(e -> e === trueminor, c.hybridedges))
	estedge = only(estchains[k].hybridedges)
	return (absterror, abs(trueminor.gamma - estedge.gamma))
end

"""
Maps a root-free key to each "chain" of `net`: a maximal path of edges whose interior nodes have
degree 2 (a root, or an old root left behind by rerooting). A chain's length is the sum of its
edge lengths, so a rooted network and its semidirected version have the same chains.

Keys ignore edge directions:
- a chain whose removal disconnects the network is keyed by the leaves on the side that does not
  contain the alphabetically first leaf
- a chain on a cycle is keyed by the pair of leaf sets hanging off its two end nodes
"""
function rootfreeedgekeys(net::HybridNetwork)
	isend(n) = length(n.edge) != 2
	other(e, n) = e.node[1] === n ? e.node[2] : e.node[1]

	chains = []
	seen = Base.IdSet{PhyloNetworks.Edge}()
	for n in net.node, e in n.edge
		(isend(n) && !(e in seen)) || continue
		cur, edge, len, hybs = n, e, 0.0, PhyloNetworks.Edge[]
		while true
			push!(seen, edge)
			len += edge.length
			edge.hybrid && push!(hybs, edge)
			cur = other(edge, cur)
			isend(cur) && break
			edge = only(x for x in cur.edge if x !== edge)
		end
		push!(chains, (a=n, b=cur, length=len, hybridedges=hybs, pendant=n.leaf || cur.leaf))
	end

	# nodes reachable from `start` without using the chains in `blocked`
	function reachable(start, blocked)
		found = Base.IdSet{PhyloNetworks.Node}([start])
		stack = [start]
		while !isempty(stack)
			n = pop!(stack)
			for (i, c) in enumerate(chains)
				(i in blocked || !(c.a === n || c.b === n)) && continue
				m = c.a === n ? c.b : c.a
				m in found || (push!(found, m); push!(stack, m))
			end
		end
		return found
	end
	reachableleaves(start, blocked) = sort([n.name for n in reachable(start, blocked) if n.leaf])

	ref = minimum(l.name for l in net.leaf)
	allleaves = sort([l.name for l in net.leaf])
	# a chain is on a cycle if its end nodes stay connected without it
	oncycle = [c.b in reachable(c.a, [i]) for (i, c) in enumerate(chains)]
	cycleidx = findall(oncycle)
	hanging(n) = reachableleaves(n, cycleidx)

	keyed = Dict{Any, Any}()
	for (i, c) in enumerate(chains)
		k = if oncycle[i]
			("cycle", sort([hanging(c.a), hanging(c.b)]))
		else
			side = reachableleaves(c.a, [i])
			("split", ref in side ? setdiff(allleaves, side) : side)
		end
		haskey(keyed, k) && error("two edges of the network share the key $k")
		keyed[k] = c
	end
	return keyed
end

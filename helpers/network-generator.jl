using PhyloNetworks, PhyloCoalSimulations, Random, StatsBase
import PhyloNetworks: addhybridedge!

function generate_network(ntaxa::Int64, nhyb::Int64, t::Float64, γ::Float64; forcetcgident::Bool=false)::HybridNetwork
	t >= 0 || error("t must be >= 0.0 (t=$t)")
	0 <= γ <= 1 || error("γ must be in [0, 1] (γ=$γ)")

	tre = readnewick(string("(", join(["t$j" for j=1:ntaxa], ","), ");"));
	for E in tre.edge E.length = 0.0 end
	tre = simulatecoalescent(tre, 1, 1)[1];

	local net::HybridNetwork
	for j = 1:1000
		net = readnewick(writenewick(tre))
		for _ = 1:nhyb
			isnothing(addhybridedge!(net, true, true; fixroot=true)) && break
		end
		net.numhybrids != nhyb && continue
		
		!forcetcgident && break
		sdnet = SNaQ.deepcopynetwork(net);
		SNaQ.semidirectnetwork!(sdnet);
		SNaQ.tcgidentifiable(sdnet) && break
	end
	net.numhybrids != nhyb && error("Could not find a valid n$(ntaxa)h$(nhyb) network after 1,000 attempts.")

	for E in net.edge E.length = t end
	for H in net.hybrid
		getparentedge(H).gamma = 1.0 - γ
		getparentedgeminor(H).gamma = γ
		getparentedgeminor(H).length = 0.0
	end
	return net
end

function parameter_errors(N::HybridNetwork, t::Float64, γ::Float64)::Vector{Float64}
	errs = Float64[];
	for E in N.edge
		getchild(E).leaf && continue
		if getchild(E).hybrid && E.ismajor
			push!(errs, E.length)
		else
			push!(errs, E.length - t)
		end
	end
	for H in N.hybrid
		γhat = min(getparentedgeminor(H).gamma, getparentedge(H).gamma)
		push!(errs, γ - γhat)
	end
	return errs
end

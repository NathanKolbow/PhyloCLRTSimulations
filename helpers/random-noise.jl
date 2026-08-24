using Distributions, Random

"""
For each row in `qCFs`:
1. Add an independent draw of `rv_distn` to each qCF
2. If the minimum qCFs (say `m`) is less than 0, add `m` to all qCFs
3. Normalize the qCFs so they sum to 1 (divide by their sum)
"""
function addqCFnoise!(qCFs::Matrix{Float64}, rv_distn::Distribution)
	for i in axes(qCFs, 1)
		qCFs[i, :] .+= rand(rv_distn, 3)
		m = min(qCFs[i, :])
		if m < 0
			qCFs[i, :] .-= m
		end
		qCFs[i, :] ./= sum(qCFs[i, :])
	end
end

@inline function runif(low::Float64, high::Float64)::Float64
	return rand() * (high - low) + low
end

@inline function runif(n::Int64, low::Float64, high::Float64)::Float64
	rvs::Vector{Float64} = Array{Float64}(undef, n)
	for i = 1:n
		rvs[i] = runif(low, high)
	end
	return rvs
end

"""
For each edge length:
1. Add an independent random draw from `rv_distn`
2. Ensure value is at least `minbl`
"""
function addedgelengthnoise!(tD::Matrix{Vector{Float64}}, rv_distn::Distribution, minbl::Float64)
	for i in axes(tD, 1)
		for j = 1:3
			for k in eachindex(tD[i, j])
				tD[i, j][k] += rand(rv_distn)
				tD[i, j][k] = max(minbl, tD[i, j][k])
			end
		end
	end
end

"""
For each row with `(n1, n2, n3)` observations for quartets 1, 2, 3:
1. Draw `n1', n2', n3'` from `Multinomial(n1+n2+n3, [n1, n2, n3])`
2. Redistribute observations according to the new `ni'` values:
	- All `ni == ni'`: skip
	- One `ni' < ni`: sample `ni' - ni` values from this quartet to put in the other columns
	- Two `ni' < ni`: put all `ni' - ni` values from these quartets into the third quartet
3. If `edgenoise_rv_distn` is not `nothing`, run `addedgelengthnoise!(tD, edgenoise_rv_distn, minbl)`
	at the end of the function

If `plusone` is true, then, `n1', n2', n3'` are drawn from `Multinomial` distribution with
probabilities `(ni + 1) / (sum(ni) + 3)` instead of `ni / sum(ni)`, thus adding a chance for
quartets with 0 observations to obtain an observation.
"""
function addjointnoise!(tD::Matrix{Vector{Float64}}, edgenoise_rv_distn::Union{Nothing, Distribution}=nothing, minbl::Float64=0.0; plusone::Bool=false)
	for i in axes(tD, 1)
		nis = [length(tD[i, 1]), length(tD[i, 2]), length(tD[i, 3])]
		ngt = sum(nis)
		newnis = !plusone ? rand(Multinomial(ngt, nis ./ ngt)) : rand(Multinomial(ngt, (nis .+ 1) ./ (ngt + 3)))
		nidiffs = newnis .- nis

		sum(nidiffs .== 0) == 3 && continue

		removedts = []
		for j in findall(nidiffs .< 0)
			ris = sort(sample(1:length(tD[i, j]), -nidiffs[j], replace=false), rev=true)
			for ri in ris
				push!(removedts, tD[i, j][ri])
				deleteat!(tD[i, j], ri)
			end
		end
		
		for j in findall(nidiffs .> 0)
			for _ = 1:nidiffs[j]
				push!(tD[i, j], removedts[1])
				deleteat!(removedts, 1)
			end
		end
	end
	!isnothing(edgenoise_rv_distn) && addedgelengthnoise!(tD, edgenoise_rv_distn, minbl)
	nothing
end

function addgtnoise!(gt::HybridNetwork, numnnis::Int64, minbl::Float64=0.0)
	for _ = 1:numnnis
		SNaQ.performrNNI1!(gt, SNaQ.samplerNNIparameters(gt, 1, Random.TaskLocalRNG())...)
	end
	
	for E in gt.edge
		E.length += rand(Normal(0, 0.01 * E.length))
		E.length = max(E.length, minbl)
	end
	return gt
end
addgtnoise!(gts::Vector{HybridNetwork}, numnnis::Int64, minbl::Float64=0.0) =
	addgtnoise!.(gts, numnnis, minbl)

anyboundaryedges(n::HybridNetwork) = any(e -> !getchild(e).leaf && !getchild(e).hybrid && e.length <= 0.01, n.edge)
anyboundarygammas(n::HybridNetwork) = any(e -> e.hybrid && e.gamma <= 0.01, n.edge)
boundaryedges(n::HybridNetwork) = [e.length for e in n.edge if !getchild(e).leaf && !getchild(e).hybrid && (e.length <= 0.01 || e.length >= 5.0)]
boundarygammas(n::HybridNetwork) = [e.gamma for e in n.edge if e.hybrid && e.gamma <= 0.01]
anyboundarygammas(n::HybridNetwork) =
	any(e -> e.hybrid && e.gamma <= 0.01, n.edge)
anyboundaries(n::HybridNetwork) = (anyboundaryedges(n), anyboundarygammas(n))
boundaries(n::HybridNetwork) = (boundaryedges(n), boundarygammas(n))

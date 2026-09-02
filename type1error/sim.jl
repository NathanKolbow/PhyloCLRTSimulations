# Simulation study to investigate the type-1 error rate (probability to reject
# H0 when H0 is true) of the CLRT and CLIC methods.
#
# Primary variables of interest:
# - Number of taxa:								8
# - Number of hybrids in the true H0 model:		0
# - Number of gene trees:						10, 100, 1000
# - edge lengths:								1.0
# - ϵ smoothing parameters:						0, 0.0001, 0.001, 0.01, 0.1, 0.5, 1.0
# Repetitions per parameter combination: 		10

Threads.nthreads() > 1 || @warn "Only $(Threads.nthreads()) threads are in use."

length(ARGS) > 0 || error("Must provide 1 command line argument for the value of ϵ (0, 0.0001, 0.001, 0.01, 0.1, 0.5, or 1.0)")
ϵ = parse(Float64, ARGS[1])
ϵ in [0, 0.0001, 0.001, 0.01, 0.1, 0.5, 1.0] || error("ϵ=$(ϵ) not allowed.")
@info "ϵ = $ϵ selected."

const NGTS  	= [10, 100, 1000];
const TS    	= [1.0];
const NREP  	= 250;
const MODELS 	= ["joint", "marginal"];
const TESTS 	= ["CLIC", "cw", "cwP", "cLR", "cLR1", "cLR2", "cLRI"];

include("../includes.jl")
outpath = joinpath(@__DIR__, "dat-eps$(ϵ).csv")
dat = isfile(outpath) ? CSV.read(outpath, DataFrame) :
	DataFrame(
		ngt=Int64[], t=Float64[], model=String[], test=String[], result=Float64[], eps=Float64[], gamma=Float64[]
	)

for irep = 1:NREP, ngt in NGTS, t in TS, model in MODELS
	hasall(dat, ϵ, model, ngt, TESTS, NREP) && continue
	print("\rmodel=$model ngt=$ngt                              ")

	# Data generation
	truenet = generate_network(8, 0, t, 0.5)	# γ doesn't matter here
	gts = simulatecoalescent(truenet, ngt, 1; inheritancecorrelation=1.0);
	dcf = gts2DCF(gts);
	tD = t_distribution(gts);
	
	# Optimizations
	print("\rmodel=$model ngt=$ngt    [r=$irep H0 SNaQ]       ")
	H0 = snaq!(readnewick(writenewick(truenet)), dcf; hmax=0, runs=100, Nfail=10, filename="")
	for E in H0.edge E.length = max(0.1, min(2.5, E.length)) end	# if H0 has boundary edges, soft-reset them
	print("\rmodel=$model ngt=$ngt    [r=$irep H1 SNaQ]       ")
	H1 = snaq!(H0, dcf; hmax=1, runs=100, Nfail=10, filename="")
	stuckid = 0
	while H1.numhybrids == 0
		stuckid += 1
		print("\rmodel=$model ngt=$ngt    [r=$irep H1 SNaQ (#$stuckid)]     ")
		H1 = snaq!(H0, dcf; hmax=1, runs=100, Nfail=10, filename="")
	end

	print("\rmodel=$model ngt=$ngt    [r=$irep H0 opt]       ")
	H0 = optimize_given_model(H0, gts, model)
	print("\rmodel=$model ngt=$ngt    [r=$irep H1 opt]       ")
	H1 = optimize_given_model(H1, gts, model)
	estγ = getparentedgeminor(H1.hybrid[1]).gamma

	# Testing
	try
		d0s, d1s = get_dists(H0, H1, model)
		lk_comps = gather_likelihood_components(d0s, d1s, gts, ϵ == 0.0 ? "rough" : "smooth", ϵ)

		for test in TESTS
			try
				push!(dat, [ngt, 1.0, model, test, run_test(H1, test, lk_comps), ϵ, estγ]; promote=true)
			catch e
			end
		end
		CSV.write(outpath, dat)
	catch
	end
end
CSV.write(outpath, dat)

Threads.nthreads() > 1 || @warn "Only $(Threads.nthreads()) threads are in use."

length(ARGS) > 0 || error("Must provide 1 command line argument for the value of ϵ (0, 0.0001, 0.001, 0.01, 0.1, or 1.0)")
ϵ = parse(Float64, ARGS[1])
ϵ in [0, 0.0001, 0.001, 0.01, 0.1, 1.0] || error("ϵ=$(ϵ) not allowed.")

const NGTS  	= [100, 10, 1000];
const TS    	= [1.0];
const GAMMAS    = 0.0:0.01:0.50;
const NREP  	= 1000;
const MODELS 	= ["joint", "marginal"];
const TESTS 	= ["cw", "cwP", "cLR", "cLR1", "cLR2", "cLRI"];

include("../includes.jl")
outpath = joinpath(@__DIR__, "dat-eps$(ϵ).csv")
dat = isfile(outpath) ? CSV.read(outpath, DataFrame) :
	DataFrame(
		simid=Int64[], gamma=Float64[], eps=Float64[], test=String[],
		model=String[], result=Float64[], t=Float64[], ngt=Int64[], 
		snaqabssumterr=Float64[], snaqabssumgammaerr=Float64[],
		newabssumterr=Float64[], newabssumgammaerr=Float64[]
	)
dat[!,:model] = convert.(String, dat[!,:model])
dat[!,:test] = convert.(String, dat[!,:test])

for irep = 1:NREP, ngt in NGTS, t in TS, model in MODELS
	simid = abs(rand(Int64))
	Random.seed!(simid)

	@simlog "\n\nngt=$ngt model=$model irep=$irep"

	truenet = generate_network(10, 1, t, 0.5; forcetcgident=true)
	for γ in GAMMAS
		@simlog "\tγ=$γ" begin
			getparentedge(truenet.hybrid[1]).gamma = 1.0 - γ
			getparentedgeminor(truenet.hybrid[1]).gamma = γ
			gts = simulatecoalescent(truenet, ngt, 1; inheritancecorrelation=1.0);

			dcf = gts2DCF(gts);
			oCFs = gts2CFs(gts);
			tD = t_distribution(gts);

			snaqH0 = snaqnetsearch(truenet, dcf, 0)
			fitnumericalparameters!(snaqH0, dcf, 1.0; maxeval=5000)
			optH0 = optimize_given_model(
				snaqH0, gts, model, ϵ; verbose=false,
				rootmaxeval=20, finalmaxeval=2500
			)
			lkH0 = gather_likelihood_components(optH0, gts, model, ϵ)
			
			snaqH1 = snaqnetsearch(truenet, dcf, 1)
			fitnumericalparameters!(snaqH1, dcf, 1.0; maxeval=5000)
			optH1 = optimize_given_model(
				snaqH1, gts, model, ϵ; verbose=false,
				rootmaxeval=20, finalmaxeval=2500
			)
			lkH1 = gather_likelihood_components(optH1, gts, model, ϵ)


			# Tests
			snaqerrors = absparamerrors(truenet, snaqH1)
			newerrors = absparamerrors(truenet, optH1)
			for test in TESTS
				push!(dat, [simid, γ, ϵ, test, model, run_test(
					optH1, test, CompTypes([lkH0[[1, 3, 4]]..., lkH1...])
				), t, ngt, snaqerrors..., newerrors...])
			end

			# CLIC
			push!(dat, [simid, γ, ϵ, "CLIC", model,
				CLICstatistic(lkH1[4], lkH1[3], lkH1[1]) - CLICstatistic(lkH0[4], lkH0[3], lkH0[1]),
				t, ngt, snaqerrors..., newerrors...
			])
		end
		CSV.write(outpath, dat)
	end
	CSV.write(outpath, dat)
end
CSV.write(outpath, dat)
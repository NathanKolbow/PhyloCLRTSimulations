@error "Try adjusting variability: sum(dlogf) * sum(dlogf)' instead of sum(dlogf * dlogf') and sensitivity sum(dlogf) * sum(dlogf)' instead of using second derivative"
Threads.nthreads() > 1 || @warn "Only $(Threads.nthreads()) threads are in use."

length(ARGS) > 0 || error("Must provide 1 command line argument for the value of ϵ (0, 0.0001, 0.001, 0.01, 0.1, 0.5, or 1.0)")
ϵ = parse(Float64, ARGS[1])
ϵ in [0, 0.0001, 0.001, 0.01, 0.1, 0.5, 1.0] || error("ϵ=$(ϵ) not allowed.")
@info "ϵ = $ϵ selected."

const NGTS  	= [10, 50, 100, 500]; #, 1000];
const TS    	= [1.0];
const NREP  	= 1000;
const MODELS 	= ["marginal", "joint"];
const TESTS 	= ["cw", "cwP", "cLR", "cLR1", "cLR2", "cLRI"];	# NOT including DDSE, Djump
const TRUEH		= [0, 1, 2] # [0, 1, 2, 3];
const HOVEREST  = 1;

SNAQONLY = length(ARGS) == 2 && parse(Bool, ARGS[2])
ϵ != 0.0 && SNAQONLY && error("To run SNaQ only set ϵ=0.0 (ϵ=$ϵ provided)")
SNAQONLY && @warn "Only running SNaQ computations and optimizations."

include("../includes.jl")
outpath = joinpath(@__DIR__, "dat-eps$(ϵ).csv")
dat = isfile(outpath) ? CSV.read(outpath, DataFrame) :
	DataFrame(
		simid=Int64[], ngt=Int64[], t=Float64[], model=String[], test=String[],
		result=Float64[], eps=Float64[], trueh=Int64[], esth=Int64[]
	)
dat[!,:model] = convert.(String, dat[!,:model])	# to stop some silly bugs w/o using `push!(...; promote=true)`
dat[!,:test] = convert.(String, dat[!,:test])	# to stop some silly bugs w/o using `push!(...; promote=true)`


for irep = 1:NREP, ngt in NGTS, t in TS, model in MODELS, trueh in TRUEH
	try
		simid = abs(rand(Int64));
		Random.seed!(simid)
		@simlog "\n\nngt=$ngt model=$model trueh=$trueh irep=$irep"

		truenet = generate_network(8, trueh, t, 0.5; forcetcgident=true)
		gts = simulatecoalescent(truenet, ngt, 1; inheritancecorrelation=1.0);
		
		dcf = gts2DCF(gts);
		oCFs = gts2CFs(gts);
		tD = t_distribution(gts);

		@simlog "\t\tInferring SNaQ networks" snaqnets = [snaqnetsearch(truenet, dcf, h) for h = 0:(trueh+HOVEREST)];
		@simlog "\t\tOptimizing networks" optnets = SNAQONLY ? [] : [
			optimize_given_model(sn, gts, model, ϵ; verbose=false, rootmaxeval=25, finalmaxeval=500)
			for sn in snaqnets
		]
		
		@simlog "\t\tGathering likelihood components" lks = SNAQONLY ? [] : [gather_likelihood_components(on, gts, model, ϵ) for on in optnets];

		@simlog "\t\tRe-fitting SNaQ networks" begin
			for sn in snaqnets
				fitnumericalparameters!(sn, dcf, 1.0; maxeval=5000)
			end
		end

		@simlog "\t\tRunning tests" begin 
			for test in TESTS
				!SNAQONLY && push!(dat, [simid, ngt, t, model, test, 0.0, ϵ, trueh, 0])
				push!(dat, [simid, ngt, t, "quartet", test, 0.0, ϵ, trueh, 0])

				for hmax = 1:(trueh+HOVEREST)
					i = hmax+1
					push!(dat, [simid, ngt, t, "quartet", test, run_test_old_model(
						snaqnets[i-1], snaqnets[i], gts, test
					), ϵ, trueh, hmax])

					SNAQONLY && continue
					push!(dat, [simid, ngt, t, model, test, run_test(
						optnets[i], test, CompTypes([lks[i-1][1], lks[i-1][3], lks[i-1][4], lks[i][1], lks[i][2], lks[i][3], lks[i][4]])
					), ϵ, trueh, hmax])
				end
			end

			# CLIC specifically
			for hmax = 0:(trueh+HOVEREST)
				i = hmax+1
				!SNAQONLY && push!(dat, [simid, ngt, t, model, "CLIC", CLICstatistic(lks[i][4], lks[i][3], lks[i][1]), ϵ, trueh, hmax])
				push!(dat, [simid, ngt, t, "quartet", "CLIC", quartetCLICstatistic(snaqnets[i], gts), ϵ, trueh, hmax])
			end
		end

		################ DDSE and Djump for new AND old models [old only if ϵ=0] ################
		if ϵ == 0.0
			@simlog "\t\tOld DDSE and Djump" try
				rows = []
				hmaxes = collect(0:(length(snaqnets)-1))
				quartet_DDSE_pens = PhyloCLRT.DDSEpenalties(SNaQscore.(snaqnets), hmaxes, ngt)
				Djumpbest = PhyloCLRT.Djumpbestmodel(SNaQscore.(snaqnets), hmaxes, ngt)
				quartet_Djump_pens = [hmax == Djumpbest ? 1.0 : 0.0 for hmax in hmaxes]

				for (DDSE_pen, Djump_pen, hmax) in zip(quartet_DDSE_pens, quartet_Djump_pens, 0:(trueh+HOVEREST))
					push!(rows, [simid, ngt, t, "quartet", "DDSE", DDSE_pen, ϵ, trueh, hmax])
					push!(rows, [simid, ngt, t, "quartet", "Djump", Djump_pen, ϵ, trueh, hmax])
				end
				# Only push rows if everything works
				for row in rows
					push!(dat, row)
				end
			catch e
			end
		end
		
		@simlog "\t\tNew DDSE and Djump" try
			SNAQONLY && error("dummy error to trigger the try-catch and exit")

			rows = []
			newLLs = SNaQscore.(optnets);
			hmaxes = collect(0:(trueh+HOVEREST))
			new_DDSE_pens = PhyloCLRT.DDSEpenalties(.-newLLs, hmaxes, ngt)
			new_Djump_best = PhyloCLRT.Djumpbestmodel(.-newLLs, hmaxes, ngt)
			new_Djump_pens = [h == new_Djump_best ? 1.0 : 0.0 for h = 0:(trueh + HOVEREST)]
			for (DDSE_pen, Djump_pen, hmax) in zip(new_DDSE_pens, new_Djump_pens, 0:(trueh+HOVEREST))
				push!(rows, [simid, ngt, t, model, "DDSE", DDSE_pen, ϵ, trueh, hmax])
				push!(rows, [simid, ngt, t, model, "Djump", Djump_pen, ϵ, trueh, hmax])
			end
			# Only push rows if everything works
			for row in rows
				push!(dat, row)
			end
		catch e
		end


		CSV.write(outpath, dat)
	catch e
		@error e
	end
end
CSV.write(outpath, dat)

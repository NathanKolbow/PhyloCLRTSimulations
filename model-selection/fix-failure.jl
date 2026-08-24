include("../includes.jl")

ϵ = 0.0
irep = 1; ngt = 50; t = 0.5; model = "joint"; trueh = 1;

optnetdict = Dict();
snaqnetdict = Dict();

const MODELS 	= ["joint", "marginal"];
const TESTS 	= ["CLIC", "cw", "cwP", "cLR", "cLR1", "cLR2", "cLRI"];	# NOT including DDSE, Djump
const TRUEH		= [0, 1, 2] #, 3];
const HOVEREST  = 1;

dat = DataFrame(
	simid=Int64[], ngt=Int64[], t=Float64[], model=String[], test=String[],
	result=Float64[], eps=Float64[], trueh=Int64[], esth=Int64[]
);


simid = abs(rand(Int64));
Random.seed!(simid);
@simlog "\n\nngt=$ngt model=$model trueh=$trueh irep=$irep"

truenet = generate_network(8, trueh, t, 0.35)
gts = simulatecoalescent(truenet, ngt, 1; inheritancecorrelation=1.0);

dcf = gts2DCF(gts);
oCFs = gts2CFs(gts);
tD = t_distribution(gts);

@simlog "\t[hmax=0]"
iterrows = []
uqhmaxes = 0
@simlog "\t\tSNaQ search:" lastnet = snaq!(majortree(truenet), dcf; hmax=0, runs=100, Nfail=10, filename="", ρ=1.0);
hypdats = [getHypothData(lastnet, oCFs, gts, 0, (a)->true)]
@simlog "\t\tTree optimization:" lastnet = optimize_given_model(lastnet, gts, model, ϵ);
snaqnet = SNaQ.deepcopynetwork(lastnet);
all_optnets = [SNaQ.deepcopynetwork(lastnet)];	# stored for DDSE and Djump
all_snaqnets = [SNaQ.deepcopynetwork(snaqnet)];

################ Data generation AND new model tests ################
# Push a test to auto-accept h=0
for test in TESTS
	push!(iterrows, [simid, ngt, t, model, test, 0.0, ϵ, trueh, 0])
end

for hmax=1:(trueh+HOVEREST)
	@simlog "\t[hmax=$hmax]"
	simpref = "ngt=$ngt model=$model trueh=$trueh irep=$irep [hmax=$hmax]"
	try
		@simlog "\t\tSNaQ search:" begin
			snaqnet = snaq!(snaqnet, dcf; hmax=hmax, runs=100, Nfail=10, filename="", ρ=1.0)
			stuckid = 0
			while snaqnet.numhybrids < hmax
				stuckid += 1
				print("\r\t\tSNaQ search [stuckid=$stuckid]:")
				snaqnet = snaq!(snaqnet, dcf; hmax=1, runs=10, Nfail=10, filename="", ρ=1.0)
				stuckid > 10_000 && error("Could not find network with snaq!.")
			end
		end

		push!(hypdats, getHypothData(snaqnet, oCFs, gts, 0, (a)->true))
		@simlog "\t\tNetwork optimization:" optnet = optimize_given_model(snaqnet, gts, model, ϵ)
		push!(all_optnets, SNaQ.deepcopynetwork(optnet));

		d0s, d1s = get_dists(lastnet, optnet, model)
		@simlog "\t\tLikelihood components:" lk_comps = gather_likelihood_components(d0s, d1s, gts, ϵ == 0.0 ? "rough" : "smooth", ϵ)
		snaqnet = optnet
		push!(all_snaqnets, snaqnet)

		uqhmaxes += 1
		@simlog "\t\tTests:" begin
			for test in TESTS
				try
					push!(iterrows, [simid, ngt, t, model, test, run_test(optnet, test, lk_comps), ϵ, trueh, hmax])
				catch e
					@error "Failed test $test"
					@error e
				end
			end
		end
	catch e
		@error e
		break
	finally
		lastnet = snaqnet
	end
end
# Only push results for tests that have the full sequence completed successfully (safeguard against unexpected errors)
nexpected = trueh + HOVEREST + 1
keep_tests = [test for test in TESTS if length([r for r in iterrows if r[5] == test]) == nexpected]
@simlog "\t\tKeeping results for $(length(keep_tests)) tests $(length(keep_tests) == length(TESTS) ? "" : "(discarding $(setdiff(TESTS, keep_tests)))")"
for r in iterrows
	r[5] in keep_tests || continue
	push!(dat, r; promote=true)
end


################ Old model tests [ONLY IF ϵ=0 SELECTED, to eliminate redundant work] ################
if ϵ == 0.0
	try
		newicks = writenewick.(gts);
		rows = []
		@simlog "\t\tOld model tests:" begin
			for hmax=1:(trueh+HOVEREST)
				testdata = TestData(hypdats[hmax], hypdats[hmax+1], -1, oCFs, newicks, ngt)
				for test in TESTS
					if hmax == 1 push!(rows, [simid, ngt, t, "quartet", test, 0.0, 0.0, trueh, 0]) end
					push!(rows, [simid, ngt, t, "quartet", test, run_test_old_model(testdata, test), 0.0, trueh, hmax])
				end
			end
		end
		# Only push rows if everything works
		for row in rows
			push!(dat, row)
		end
	catch e
	end
end

################ DDSE and Djump for new AND old models [old only if ϵ=0] ################
if ϵ == 0.0
	@simlog "\t\tOld DDSE and Djump" try
		rows = []
		old_DDSE_pens = DDSEtest(hypdats, ngt)[2]
		old_Djump_pens = Djumptest(hypdats, ngt)[2]

		for (DDSE_pen, Djump_pen, hmax) in zip(old_DDSE_pens, old_Djump_pens, 0:(trueh+HOVEREST))
			push!(rows, [simid, ngt, t, "quartet", "DDSE", DDSE_pen, 0.0, trueh, hmax])
			push!(rows, [simid, ngt, t, "quartet", "Djump", Djump_pen, 0.0, trueh, hmax])
		end
		# Only push rows if everything works
		for row in rows
			push!(dat, row)
		end
	catch e
	end
end

@simlog "\t\tNew DDSE and Djump" try
	rows = []
	newLLs = SNaQscore.(all_optnets);
	new_DDSE_pens = DDSEpenalties(.-newLLs, collect(0:(trueh+HOVEREST)), ngt)
	new_Djump_best = Djumpbestmodel(.-newLLs, collect(0:(trueh+HOVEREST)), ngt)
	new_Djump_pens = abs.(1.0 .- [h == new_Djump_best ? 1.0 : 0.0 for h = 0:(trueh + HOVEREST)])
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

optnetdict[ngt] = SNaQ.deepcopynetwork.(all_optnets);
snaqnetdict[ngt] = SNaQ.deepcopynetwork.(all_snaqnets);

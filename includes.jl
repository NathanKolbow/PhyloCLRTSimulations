# Includes all relevant CLRT source files
B = @__DIR__;

@info "Gathering packages and include files. Output hidden until complete."
so = stdout; se = stderr;
redirect_stdout(devnull)
redirect_stderr(devnull)

using Pkg
Pkg.activate(B);
using DataFrames, CSV, RCall	# RCall needs to be loaded before the block containing `@rlibrary` for some reason...
using PhyloNetworks, SNaQ, InPhyNet, PhyloCoalSimulations, PhyloCLRT
import PhyloCLRT: TestData, HypothesisData, optimize_parameters∇!, optimize_root_placement,
	getHypothData,
	mdm_optimize_parameters∇!, mdm_optimize_root_placement,
	logf, ∇logf, sensitivity, variability,
	smooth_logf, smooth_∇logf, smooth_sensitivity, smooth_variability,
	cw, cwP, CLICstatistic, cLR, cLR1, cLR2, cLRI

using LinearAlgebra, StatsBase, Random

try
	# incl_files = [
	# 	joinpath(B, "../CLRT/src/model_expansion/hypothesis-tests.jl"),
	# 	joinpath(B, "../CLRT/src/value_computations.jl"),
	# 	joinpath(B, "../CLRT/data/data_structs.jl"),
	# 	joinpath(B, "../CLRT/data/data_generation.jl"),
	# 	joinpath(B, "../CLRT/hypothesis-tests/capushe.jl"),
	# 	joinpath(B, "../CLRT/hypothesis-tests/cLRI.jl"),
	# 	joinpath(B, "../CLRT/hypothesis-tests/cwP.jl"),
	# 	joinpath(B, "../CLRT/hypothesis-tests/clic.jl")
	# ]
	# for file in incl_files
	# 	include(file)
	# end

	incl_folders = [joinpath(B, "helpers")];
	for dir in incl_folders
		for file in readdir(dir; join=true)
			include(file)
		end
	end

	# Load the libs while output is hidden
	R"""
	library(capushe)
	library(tidyverse)
	"""
finally
	redirect_stdout(so)
	redirect_stderr(se)
end
@info "Completed."
println()

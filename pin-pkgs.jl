# Locks package versions in Manifest.toml
# and sets the [compat] section of the
# Project.toml so that current package
# versions are preserved for reproducibility.

using Pkg
Pkg.activate(@__DIR__)
proj = Pkg.project()
deps = proj.dependencies
topin = [Pkg.PackageSpec(k, v) for (k, v) in deps]

for pkg in topin
	Pkg.pin(pkg)
end


using Pkg, TOML
direct_deps = [
    pkg for pkg in values(Pkg.dependencies()) 
    if pkg.is_direct_dep && pkg.version !== nothing
]
project_file = Base.active_project()
project_toml = TOML.parsefile(project_file)
if !haskey(project_toml, "compat")
    project_toml["compat"] = Dict{String, Any}()
end
for dep in direct_deps
    project_toml["compat"][dep.name] = "=$(dep.version)"
end
open(project_file, "w") do io
    TOML.print(io, project_toml, sorted=true)
end

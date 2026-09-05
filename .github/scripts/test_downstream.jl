using Pkg
using TOML

qi_path = ENV["QI_PATH"]
downstream_path = ENV["DOWNSTREAM_PATH"]
package = TOML.parsefile(joinpath(downstream_path, "Project.toml"))["name"]

Pkg.develop(path = qi_path)
Pkg.develop(path = downstream_path)
Pkg.build(package)
Pkg.test(package)

using Pkg

Pkg.develop(path = ENV["QI_PATH"])
for path in readdir(ENV["DOWNSTREAMS_PATH"]; join = true)
    Pkg.develop(path = path)
end

modules = map((
    :QuantumInterface,
    :QuantumOpticsBase,
    :QuantumOptics,
    :SecondQuantizedAlgebra,
    :Gabs,
    :QuantumSymbolics,
    :QuantumClifford,
    :QuantumSavory,
)) do name
    Base.eval(Main, :(import $name))
    getfield(Main, name)
end
import TensorCore

for mod in modules
    for name in (:tensor, :⊗)
        isdefined(mod, name) || continue
        @assert getfield(mod, name) === TensorCore.tensor "$mod.$name is not TensorCore.tensor"
    end
end
@assert !applicable(TensorCore.tensor)

ambiguities = Base.detect_ambiguities(modules..., TensorCore; recursive = true)
tensor_ambiguities = filter(ambiguities) do (left, right)
    any((left, right)) do method
        isdefined(method.module, method.name) &&
            getfield(method.module, method.name) === TensorCore.tensor
    end
end
@assert isempty(tensor_ambiguities) "TensorCore.tensor ambiguities: $tensor_ambiguities"

module GridapMT

using Gridap
import Gridap: solve
using Gridap.Fields

using GridapGmsh
import GridapGmsh: gmsh

using Gridap.Geometry
using Gridap.CellData
using Gridap.Visualization
using Gridap.FESpaces

using JLD2

abstract type ForwardProblem end

mutable struct MT <: ForwardProblem
    model::UnstructuredDiscreteModel
    freqs::Vector{Float64}
    obs::Vector{Float64}
    # σ₀::Float64
end

mutable struct FEProblem
    order::Int
    Vₑ
    Vₕ
    V₀
    Ω
    dΩ
end

mutable struct ParameterFields
    σ::FEFunction
    ρ::FEFunction
end

export MT, FEProblem, ParameterFields

include("dirichletBC.jl")
export dirichlet_hom_E, dirichlet_hom_H

include("gmsh.jl")
export loadModel, PrismGenerator, gmsh

include("rhoa.jl")
export get_rhoa_phase

include("utils.jl")
export getPT

include("fem.jl")
export FE_setup

include("solve.jl")
export solve_eh, solve_e, solve_h

include("grid.jl")
export TensorGrid
export createTensorGrid, createGrid, interpolate_grid_to_fe_space

end # module GridapMT

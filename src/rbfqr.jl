module rbfqr

using LinearAlgebra

struct Psistruct
	j::Vector{Int64}
	m::Vector{Int64}
  p::Vector{Int64}
	cs::Vector{Int64}
  ep::Float64
	xk::Matrix{Float64}        
  Rt::Matrix{Float64}
	columns::Vector{Float64}
	rr::Float64
	cc::Matrix{Float64}
	A0::Matrix{Float64}
end

include("Discretization.jl")
include("rbfqr_diffmat_2d.jl")
include("init_psi_2d.jl")

export Discretization
export rbfqr_diffmat_2d
export Psistruct, Tstruct, Qstruct, Rstruct


end # end



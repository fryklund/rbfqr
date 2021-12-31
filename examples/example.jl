using rbfqr
using LinearAlgebra
using rbfqr.RBFdiscretization

include("../tests/testfunction.jl")

# The shape parameter
ep = 0.5

# Use Halton nodes in the unit square
N = 200
d=2

xk = halton(N,d)

# Evaluation points
Ne = 40 # Per dimension
x = range(0,1,length=Ne)
xx,yy = ndgrid(x,x)
xe = hcat(xx[:], yy[:])
A, Psi = rbfqr_diffmat_2d("1", xe, xk, ep)
Ax = rbfqr_diffmat_2d("x", xe, Psi)[1]
Ay = rbfqr_diffmat_2d("y", xe, Psi)[1]
Axx = rbfqr_diffmat_2d("xx", xe, Psi)[1]
Axy = rbfqr_diffmat_2d("xy", xe, Psi)[1]
Ayy = rbfqr_diffmat_2d("yy", xe, Psi)[1]
L = rbfqr_diffmat_2d("L", xe, Psi)[1]

uk = f("0",xk)

diff = A*uk-f("0",xe)
diffx = Ax*uk-f("x",xe)
diffy = Ay*uk-f("y",xe)
diffxx = Axx*uk-f("xx",xe)
diffxy = Axy*uk-f("xy",xe)
diffyy = Ayy*uk-f("yy",xe)
diffL = L*uk-f("L",xe)

println("max error: ", maximum(diff))
println("max error x: ", maximum(abs.(diffx)))
println("max error y: ", maximum(abs.(diffy)))
println("max error xx: ", maximum(abs.(diffxx)))
println("max error xy: ", maximum(abs.(diffxy)))
println("max error yy: ", maximum(abs.(diffyy)))
println("max error L: ", maximum(abs.(diffL)))

#=
# Not yet implemented
op = ["1", "x", "y", "xx", "xy", "yy", "L"]
Aall, Psi = rbfqr_diffmat_2d(op, xe, xk, ep); 
=#

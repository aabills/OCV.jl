using Test
using OCV

using Test
using OCV
using BenchmarkTools

A = [-123466.90208905171,
-65215.58255076725,
-77739.69488015345,
-132750.8136541972,
-109183.62667589705]

U_0 = 3.4271972387993173      # U_0_p   # U_0_n

c_s_max = 100000.0
c_s_min = 0.0

cathodeocv = OCV.RKPolynomial(A,U_0,c_s_max,c_s_min,length(A_p))

atol = 1e-10

x = 0.56
V = cathodeocv(x,285.0)

x̂ = OCV.mybisection(V,0.6,0.4,cathodeocv,285.0,atol=atol)
V̂ = cathodeocv(x̂,285.0)

@test x ≈ x̂
@test V ≈ V̂




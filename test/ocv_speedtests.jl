using Test
using OCV
using BenchmarkTools

A_p = [-123466.90208905171,
-65215.58255076725,
-77739.69488015345,
-132750.8136541972,
-109183.62667589705]

A_n = [-9.04729951400691e6,
-8.911003506471587e6,
-9.04657963355355e6,
-8.904669509837592e6,
-9.0363250622869e6,
-9.05878665345401e6,
-9.606335422964232e6,
-8.023042975317075e6,
-2.3190474522951595e6,
 1.4303914546788693e6]

U_0_p = 3.4271972387993173      # U_0_p
U_0_n = -46.09780594535385    # U_0_n

c_s_max⁺ = 100000.0
c_s_max⁻ = 100000.0
c_s_min⁺ = 0.0
c_s_min⁻ = 0.0

cathodeocv = OCV.RKPolynomial(A_p,similar(A_p),U_0_p,c_s_max⁺,c_s_min⁺,length(A_p))
anodeocv   = OCV.RKPolynomial(A_n,similar(A_n),U_0_n,c_s_max⁻,c_s_min⁻,length(A_n))

x = 0.2
benchmark = @benchmark OCV.calcocv($anodeocv,$x,285.0)
V = OCV.calcocv(anodeocv,x,285.0)


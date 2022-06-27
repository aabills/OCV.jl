module OCV
using JuMP
import Ipopt

#
# Type to Represent an RK Polynomial
#

abstract type OpenCircuitVoltage end

#Evaluate OpenCircuitVoltage types
(c::OpenCircuitVoltage)(x,T) = calcocv(c,x,T)
(c::OpenCircuitVoltage)(x) = calcocv(c,x,297.)

struct RKPolynomial{T,S} <: OpenCircuitVoltage
    Aterms::T
    U_0::S
    c_s_max::S
    nA::Int
    # UpperFillingFraction
    # LowerFillingFraction
end

#
# Constructor for RK Polynomial
#
function RKPolynomial(Aterms,U_0)
    return RKPolynomial(Aterms,U_0,length(Aterms))
end

#
# Function to Calculate RK Polynomials
#
function calcocv(RK::RKPolynomial,x::W,T) where {W}
    F=96485.33212
    n=1.
    R=8.314
    a=1
    rk::W = RK.Aterms[1]*((2 *x-1)^((a-1) +1))
    # Activity Correction Summation
    @inbounds for a in 2:RK.nA
       rk += RK.Aterms[a]*((2 *x-1)^((a-1) +1) - (2 *x*(a-1) *(1 -x)/ (2 *x-1)^(1 -(a-1))))
    end
    voltage::W = rk/(n*F)+RK.U_0+((R*T)/(n*F) * log((1 -x)/x))
    return voltage
end


# 
# Functions to force monotonicity in OpenCircuitVoltage Fits
# 
function MonotonicIncreaseLeastSquaresFit(A,y)
    model = Model(Ipopt.Optimizer)
    dim = length(y)
    @variable(model,x[1:dim])
    @variable(model,ŷ[1:dim])
    for i in 2:dim
        @constraint(model,ŷ[i]>=ŷ[i-1])
    end
    @constraint(model,ŷ.==A*x)
    @objective(model,Min,(ŷ.-y)'*((ŷ.-y)))
    optimize!(model)
    return model,x,ŷ
end
# 
# Functions to force monotonicity in OpenCircuitVoltage Fits
# 
function MonotonicDecreaseLeastSquaresFit(A,y)
    model = Model(Ipopt.Optimizer)
    dim = length(y)
    @variable(model,x[1:dim])
    @variable(model,ŷ[1:dim])
    for i in 2:dim
        @constraint(model,ŷ[i]<=ŷ[i-1])
    end
    @constraint(model,ŷ.==A*x)
    @objective(model,Min,(ŷ.-y)'*((ŷ.-y)))
    optimize!(model)
    return model,x,ŷ
end

export RKPolynomial,calcocv,MonotonicIncreaseLeastSquaresFit,MonotonicDecreaseLeastSquaresFit

end # module

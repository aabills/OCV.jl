module OCV
using Reexport
@reexport using JuMP
@reexport using DataInterpolations
import Ipopt

const F=96485.33212
const n=1.
const R=8.314



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
    c_s_min::S
    nA::Int
    # UpperFillingFraction
    # LowerFillingFraction
end

#
# Constructor for RK Polynomial
#
function RKPolynomial(Aterms,U_0,c_s_max,c_s_min)
    return RKPolynomial(Aterms,U_0,c_s_max,c_s_min,length(Aterms))
end

#
# Function to Calculate RK Polynomials
#
function calcocv(RK::RKPolynomial,x::W,T) where {W}
    F=96485.33212
    n=1.
    R=8.314
    a=1
    @inbounds rk::W = RK.Aterms[1]*((2 *x-1)^((a-1) +1))
    # Activity Correction Summation
    @simd for a in 2:RK.nA
       @inbounds @fastmath rk += RK.Aterms[a]*((2 *x-1)^((a-1) +1) - (2 *x*(a-1) *(1 -x)/ (2 *x-1)^(1 -(a-1))))
    end
    voltage::W = @fastmath rk/(n*F)+RK.U_0+((R*T)/(n*F) * log((1 -x)/x))
    return voltage
end

#
# Function to calculate the filling fraction of an OCV given a voltage
# 
function get_x_from_voltage(RK::OpenCircuitVoltage, V, T;atol=1e-10)
# Find if this is an increasing or decreasing OCV
x₁ = 0.2
x₂ = 0.3
V₁ = calcocv(RK,x₁,T)
V₂ = calcocv(RK,x₂,T)
if V₁ ≥ V₂
    return OCV.mybisection(V,1.0,0.0,RK,T,atol=atol)
else
    return OCV.mybisection(V,0.0,0.1,RK,T,atol=atol)
end

end

function RK_Matrix(x, num_terms)
    mat = zeros(length(x), num_terms+1)
    a=1
    mat[:,1] .= ((2 .*x.-1).^((a.-1) .+1))
    # Activity Correction Summation
    for a in 2:num_terms
         mat[:, a] = ((2 .*x.-1).^((a.-1) .+1) .- (2 .*x.*(a.-1) .*(1 .-x)./ (2 .*x.-1).^(1 .-(a.-1))))
    end
    mat[:, end] .= OCV.F
    return mat
end

function get_lhs(voltage, x, T)
    if length(voltage) != length(x)
        error("length of voltage must equal length of x")
    end
    lhs = n*OCV.F.*voltage - OCV.R.*T.*log.((1 .-x)./x)
    return lhs
end





#
# Quick and dirty bisection implementation
#
function mybisection(V, min, max, RK::OpenCircuitVoltage, T; atol=1e-10)
    mid = (min + max) / 2
    V̂ = calcocv(RK,mid,T)
    if abs(V̂-V)≤atol
        return mid
    elseif V̂ > V
        return mybisection(V,min,mid,RK,T,atol=atol)
    else
        return mybisection(V,mid,max,RK,T,atol=atol)
    end
end


# 
# Functions to force monotonicity in OpenCircuitVoltage Fits
# 
function MonotonicIncreaseLeastSquaresFit(A,y)
    model = Model(Ipopt.Optimizer)
    num_As = size(A)[2]
    dim = length(y)
    @variable(model,x[1:num_As])
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
    num_As = size(A)[2]
    dim = length(y)
    @variable(model,x[1:num_As])
    @variable(model,ŷ[1:dim])
    for i in 2:dim
        @constraint(model,ŷ[i]<=ŷ[i-1])
    end
    @constraint(model,ŷ.==A*x)
    @objective(model,Min,(ŷ.-y)'*((ŷ.-y)))
    optimize!(model)
    return model,x,ŷ
end

export RKPolynomial,calcocv,MonotonicIncreaseLeastSquaresFit,MonotonicDecreaseLeastSquaresFit,get_x_from_voltage

include("spline_ocvs.jl")
export SplineOCV


end # module

struct SplineOCV{W} <: OpenCircuitVoltage
    spline
    c_s_max::W
    c_s_min::W
end

#
# Spline Calculate OCV. Currently no support for Temperature changes
#
function calcocv(S::SplineOCV, x::W, T) where {W}
    return S.spline(x)
end

#
# Constructor Function
#
function SplineOCV(x, V, splineType, c_s_max, c_s_min; kwargs...)
    spline = splineType(V,x,kwargs...)
    return SplineOCV(spline, c_s_max, c_s_min)
end





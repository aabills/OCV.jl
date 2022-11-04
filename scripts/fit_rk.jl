using OCV


function RK_Matrix(xs,nc)
    F=96485.33212        #Faraday's Constant [Coulombs/mol]
    n=1 
    R=8.314
    nr = length(xs)
    RK_mat = zeros(nr,nc)
    a=1
    #First Col is different
    @. RK_mat[:,1] = ((2 *xs-1)^((a-1) +1))
    # Activity Correction Summation
    for a in 2:nc
       @. RK_mat[:,a] = ((2 *xs-1)^((a-1) +1) - (2 *xs*(a-1) *(1 -xs)/ (2 *xs-1)^(1 -(a-1))))
    end
    RK_mat = RK_mat./(n.*F)
    RK_mat = hcat(RK_mat,ones(nr))
    return RK_mat
end

xs = collect(range(0.1,stop=0.9,length=1000))
ys = collect(1:1000)


A = RK_Matrix(xs,5)

a = OCV.MonotonicIncreaseLeastSquaresFit(A,ys)


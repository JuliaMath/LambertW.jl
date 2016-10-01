export lamwcoeff, cmpbp, cmpfp, xlambertw, bigbp
export cmpxlam,cmpexp

# Testing and development for expansion about branch point.  Refer to
# the paper "On the Lambert W function".  In (4.22) coefficients μ₀
# through μ₃ are given explicitly. Recursion relations (4.23) and
# (4.24) for all μ are also given. This code implements the recursion
# relations.

# (4.23) and (4.24) give zero based coefficients
cset(a,i,v) = a[i+1] = v
cget(a,i) = a[i+1]

# (4.24)
function compa(k,m,a)
    sum0 = zero(eltype(m))
    for j in 2:k-1
        sum0 += cget(m,j) * cget(m,k+1-j)
    end
    cset(a,k,sum0)
    sum0
end

# (4.23)
function compm(k,m,a)
    kt = convert(eltype(m),k)
    mk = (kt-1)/(kt+1) *(cget(m,k-2)/2 + cget(a,k-2)/4) -
        cget(a,k)/2 - cget(m,k-1)/(kt+1)
    cset(m,k,mk)
    mk
end

# We plug the known value μ₂ == -1//3 for (4.22) into (4.23) and
# solve for α₂. We get α₂ = 0.
# compute array of coefficients μ in (4.22).
# m[1] is μ₀
function lamwcoeff(T::DataType, n::Int)
    a = Array(T,n)
    m = Array(T,n)
    cset(a,0,2)  # α₀ literal in paper
    cset(a,1,-1) # α₁ literal in paper
    cset(a,2,0)  # α₂ get this by solving (4.23) for alpha_2 with values printed in paper
    cset(m,0,-1) # μ₀ literal in paper
    cset(m,1,1)  # μ₁ literal in paper
    cset(m,2,-1//3) # μ₂ literal in paper, but only in (4.22)
#    compa(3,m,a)  # α₃, from (4.24) get 1/9,   who knows if correct
#    compm(3,m,a)  # μ₃, from (4.23) get 11 / 72  agrees with explicit number in paper
#                  and this relies on the value α₃
    for i in 3:n-1  # coeffs are zero indexed
        compa(i,m,a)
        compm(i,m,a)
    end
    return m
end

const LAMWMU_FLOAT64 = lamwcoeff(Float64,100)
const LAMWMU_BIGRAT = lamwcoeff(Rational{BigInt},20)
const LAMWMU_BIGFLOAT = lamwcoeff(BigFloat,100)
const LAMWMU_RAT = lamwcoeff(Rational{Int},10)

function wser{T<:AbstractArray}(p,a::T,n::Int)
    sum0 = zero(p)
    for i in 2:n
        sum0 += p^(i-1) * a[i]
    end
    sum0
end


function _lambertw0(x::Float64,n) # 1 + W(-1/e + x)  , k = 0
    ps = 2*convert(typeof(x),e)*x;
    p = sqrt(ps)
    wser(p,LAMWMU_FLOAT64,n)
end

function _lambertwm1(x::Float64,n) # 1 + W(-1/e + x)  , k = -1
    ps = 2*convert(typeof(x),e)*x;
    p = -sqrt(ps)
    wser(p,LAMWMU_FLOAT64,n)
end

function _lambertw0(x::BigFloat,n) # 1 + W(-1/e + x)  , k = 0
    ps = 2*e*x;
    p = sqrt(ps)
    wser(p,LAMWMU_BIGFLOAT,n)
end

function _lambertwm1(x::BigFloat,n) # 1 + W(-1/e + x)  , k = -1
    ps = 2*e*x;
    p = -sqrt(ps)
    wser(p,LAMWMU_BIGFLOAT,n)
end

function alambertwbp{T<:Real}(x::T,k::Int,n::Int)
    k == 0 && return _lambertw0(x,n)
    k == -1 && return _lambertwm1(x,n)
    error("exansion about branch point only implemented for k = 0 and -1")
end

# do a big float root finding computation equivalent to alambertwbp
function bigbp(x,k)
    one(x) + lambertw(-one(x)/convert(typeof(x),e) + x, k)
end

# compare error from BigFloat exact (root finding) and expansion about branch point
# x = distance from branch point
# k = branch index, either 0 or -1
# n = order of expansion (n-1 terms because we omit the constant -1 )
# The result is that n > 5 no longer changes the result in the best case.
# Note that radius of convergence is reported to be sqrt(2), we are well within this.
# So, it makes sense to include these terms, but no more:
#    μ₄ = -113//1080
#    μ₅ =  697//17280
#
#    For smaller x, fewer coefficients are needed.
function cmpbp(x,k,n)
    w1 = alambertwbp(x,k,n)
    w2 = bigbp(x,k)
    float64(w1-w2)
end

xlambertw(x,k) = one(x) + lambertw(-one(x)/e + 1/x,k)

# The following three functions are good for comparing accuracy
# of expansion vs root finding near the bp.
# This confirms that the crossover is still about about 1/x = 10^7
cmpxlam(x,k) = (a=xlambertw(BigInt(x),k); (xlambertw(x,k) - a)/a)
cmpexp(x,k) = (a=xlambertw(BigInt(x),k); (alambertwbp(1/float(x),k) -a)/a)
cmpexp(x,k,n) = (a=xlambertw(BigInt(x),k); (alambertwbp(1/float(x),k,n) -a)/a)

function cmpfp(x,k)
    w1 = lambertw(-1/e + x,k)
    w2 = lambertw(-1/big(e)+x,k)
    float64(w1-w2)
end

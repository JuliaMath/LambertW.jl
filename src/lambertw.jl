export lambertw, lambertwbp
export ω, omega

import Base: convert

## Lambert W function

# Use Halley's iterative method to find x = lambertw(z)
# with initial point x.
function _lambertw{T<:Number}(z::T, x::T)
    two_t = convert(T,2)    
    lastx = x
    lastdiff = 0.0
    for i in 1:100
        ex = exp(x)
        xexz = x * ex - z
        x1 = x + 1
        x = x - xexz / (ex * x1 - (x + two_t) * xexz / (two_t * x1 ) )
        diff = abs(lastx - x)
        diff <= eps(abs(lastx)) && break
        if lastdiff == diff
            warn("lambertw did not converge. diff=$diff")
            break
        end
        lastx = x
        lastdiff = diff
    end
    x
end

# Real x, k = 0
# fancy initial condition does not seem to help speed.
function lambertwk0{T<:Real}(x::T)
    one_t = one(T)
    two_t = convert(T,2)
    x < -one_t/convert(T,e) && return NaN
    if x > one_t
        lx = log(x)
        llx = log(lx)
        x1 = lx - llx - log(one_t - llx/lx)/two_t
    else
        x1 = 0.567 * x
    end
    _lambertw(x,x1)
end

# Real x, k = -1
function _lambertwkm1{T<:Real}(x::T)
    -one(T)/convert(T,e) < x < zero(T) || return NaN
    _lambertw(x,log(-x))
end

function lambertw(z::Complex, k::Int)
    rT = typeof(real(z))
    one_t = one(rT)    
    if abs(z) <= one_t/convert(rT,e)
        if z == 0
            k == 0 && return z
            return -inf(rT)
        end
        if k == 0
            w = z
        elseif k == -1 && imag(z) == 0 && real(z) < 0
            w = log(-real(z))
        else
            w = log(z)
            k != 0 ? w += complex(0,k * 2 * pi) : nothing            
        end
    elseif k == 0 && imag(z) <= 0.7 && abs(z) <= 0.7
        w = abs(z+0.5) < 0.1 ? imag(z) > 0 ? complex(0.7,0.7) : complex(0.7,-0.7) : z
    else
        if real(z) == inf(rT)
            k == 0 && return z
            return z + complex(0,2*k*pi)
        end
        real(z) == -inf(rT) && return -z + complex(0,(2*k+1)*pi)
        w = log(z)
        k != 0 ? w += complex(0, 2*k*pi) : nothing
    end
    _lambertw(z,w)
end

function lambertw{T<:Real}(x::T, k::Int)
    k == 0 && return lambertwk0(x)
    k == -1 && return _lambertwkm1(x)
    error("lambertw: real x must have k == 0 or k == -1")
end

lambertw(z::Complex{Int}, k::Int) = lambertw(float(z),k)

function lambertw{T<:Integer}(x::T, k::Int)
    if k == 0
        x == 0 && return zero(x)
        x == 1 && return SpecFun.omega
    end
    lambertw(float(x),k)
end

lambertw(::MathConst{:e}) = 1
function lambertw(::MathConst{:e}, k::Int)
    k == 0 && return 1
    k == -1 && return NaN
    error("lambertw: real x must have k == 0 or k == -1")    
end
    
lambertw(x::Number) = lambertw(x,0)

## omega constant

# These literals have more than Float64 and BigFloat 256 precision
const omega_const_ = 0.567143290409783872999968662210355
const omega_const_bf_ = BigFloat("0.5671432904097838729999686622103555497538157871865125081351310792230457930866845666932194")

# maybe compute higher precision. converges very quickly
function omega_const(::Type{BigFloat})
    get_bigfloat_precision() <= 256 && return omega_const_bf_
    myeps = eps(BigFloat)
    oc = omega_const_bf_
    for i in 1:100
        nextoc = (1 + oc) / (1 + exp(oc))
        abs(oc - nextoc) <= myeps && break
        oc = nextoc
    end
    oc
end

const ω = MathConst{:ω}()
const omega = ω
convert(::Type{BigFloat}, ::MathConst{:ω}) = omega_const(BigFloat)
convert(::Type{Float64}, ::MathConst{:ω}) = omega_const_


## Expansion about branch point x = -1/e

# Better to compute only necessary terms, but this
# requires some logic.
wser(p,ps) = p - ps / 3 + (11/72) * p * ps

function _lambertw0(x) # 1 + W(-1/e + x)  , k = 0
    ps = 2*e*x;
    p = sqrt(ps)
    wser(p,ps)
end

function _lambertwm1(x) # 1 + W(-1/e + x)  , k = -1
    ps = 2*e*x;
    p = -sqrt(ps)
    wser(p,ps)
end

function lambertwbp{T<:Real}(x::T,k::Int)
    k == 0 && return _lambertw0(x)
    k == -1 && return _lambertwm1(x)
    error("exansion about branch point only implemented for k = 0 and -1")
end

lambertwbp(x) = lambertwbp(x,0)

# note that pz(z) uses the difference between numbers that are nearly equal.
# This is suggested in the paper "On the Lambert W function". But, it is not a good idea.
# Instead we use pzdiff above
# Series near branch point z = -1/e
#pz(z) = -sqrt(2*(e*z+1))
# wexp is series expansion about -1/e,  1 + W(-1/e + x)
#wexp(x) = wser(pz(x))
# wexpdiff(z) computes lambertw(-1/e+x,-1) for small positive x

### Following is only playing with rewriting expressions.

iscall(ex::Expr) = ex.head == :call
isop(ex::Expr, s::Symbol) = iscall(ex) && ex.args[1] == s
function is_x_times_fofx (ex::Expr, f::Symbol)
    isop(ex, :*) || return false
    a = ex.args
    length(a) == 3 || return false
    t = a[3]
    isop(t, f) || return false
    x = t.args[2]
    x == a[2] && return x
    return false
end

# This is an interesting toy. But there is no mechanism to
# reduce expressions to normal forms.
# So we don't know that exp(-1) == 1/e, etc.
function lambertw(ex::Expr, k::Int)
    ex == :(-pi/2) && return : (complex(0,pi/2))
    ex == :(-1/e) && return -1
    ex == :(exp(-1)) && return -1
    res = is_x_times_fofx(ex,:exp)
    res != false && return res
    res = is_x_times_fofx(ex,:log)
    res != false && return :(log($res))
    :(lambertw($ex,$k))
end

lambertw(ex::Expr) = lambertw(ex,0)

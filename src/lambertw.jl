export lambertw, lambertwbp
export ω, omega

import Base: convert

#### Lambert W function ####

# Maybe finish implementing this later ?
lambert_verbose() = false 

# Use Halley's root-finding method to find x = lambertw(z) with
# initial point x.
function _lambertw{T<:Number}(z::T, x::T)
    two_t = convert(T,2)    
    lastx = x
    lastdiff = 0.0
    for i in 1:100
        ex = exp(x)
        xexz = x * ex - z
        x1 = x + 1
        x = x - xexz / (ex * x1 - (x + two_t) * xexz / (two_t * x1 ) )
#        x = x - 2*x1*xexz/(2*x1*x1*ex-xexz*(x1+two_t))  slower than line above
        xdiff = abs(lastx - x)
        xdiff <= 2*eps(abs(lastx)) && break
        if lastdiff == diff
            lambert_verbose() && warn("lambertw did not converge. diff=$xdiff")
            break
        end
        lastx = x
        lastdiff = xdiff
    end
    x
end

### Real z ###

# Real x, k = 0
# fancy initial condition does not seem to help speed.
function lambertwk0{T<:Real}(x::T)
    const one_t = one(T)    
    const oneoe = -one_t/convert(T,e)
    x == oneoe && return -one_t
    const itwo_t = 1/convert(T,2)    
    oneoe <= x || throw(DomainError())
    if x > one_t
        lx = log(x)
        llx = log(lx)
        x1 = lx - llx - log(one_t - llx/lx) * itwo_t
    else
        x1 = 0.567 * x
    end
    _lambertw(x,x1)
end

# Real x, k = -1
function _lambertwkm1{T<:Real}(x::T)
    const oneoe = -one(T)/convert(T,e)
    x == oneoe && return -one(T)
    oneoe <= x || throw(DomainError())
    x == zero(T) && return -convert(T,Inf)
    x < zero(T) || throw(DomainError())
    _lambertw(x,log(-x))
end

function lambertw{T<:Real}(x::T, k::Int)
    k == 0 && return lambertwk0(x)
    k == -1 && return _lambertwkm1(x)
    throw(DomainError())  # more informative message like below ?
#    error("lambertw: real x must have k == 0 or k == -1")
end

function lambertw{T<:Integer}(x::T, k::Int)
    if k == 0
        x == 0 && return float(zero(x))
        x == 1 && return convert(typeof(float(x)),SpecFun.omega) # must be more efficient way
    end
    lambertw(float(x),k)
end

### Complex z ###

# choose initial value inside correct branch for root finding
function lambertw(z::Complex, k::Int)
    rT = typeof(real(z))
    one_t = one(rT)    
    if abs(z) <= one_t/convert(rT,e)
        if z == 0
            k == 0 && return z
            return complex(-convert(rT,Inf),zero(rT))
        end
        if k == 0
            w = z
        elseif k == -1 && imag(z) == 0 && real(z) < 0
            w = complex(log(-real(z)),1/10^7) # need offset for z ≈ -1/e.
        else
            w = log(z)
            k != 0 ? w += complex(0,k * 2 * pi) : nothing            
        end
    elseif k == 0 && imag(z) <= 0.7 && abs(z) <= 0.7
        w = abs(z+0.5) < 0.1 ? imag(z) > 0 ? complex(0.7,0.7) : complex(0.7,-0.7) : z
    else
        if real(z) == convert(rT,Inf)
            k == 0 && return z
            return z + complex(0,2*k*pi)
        end
        real(z) == -convert(rT,Inf) && return -z + complex(0,(2*k+1)*pi)
        w = log(z)
        k != 0 ? w += complex(0, 2*k*pi) : nothing
    end
    _lambertw(z,w)
end

lambertw(z::Complex{Int}, k::Int) = lambertw(float(z),k)

# lambertw(e + 0im,k) is ok for all k
function lambertw(::MathConst{:e}, k::Int)
    k == 0 && return 1
    throw(DomainError())
end

lambertw(x::Number) = lambertw(x,0)

### omega constant ###

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


### Expansion about branch point x = -1/e  ###

# Better to omit terms that don't contribute significant digits. But,
# that requires some logic. We get ps for free because it is needed to
# compute p. The entire call to lambertwbp(x,k) for three terms takes
# about 6 ns on my machine.
 #        -1//1          
 #         1//1          
 #        -1//3          
 #        11//72         
 #       -43//540        
 #       769//17280      
 #      -221//8505       
 #    680863//43545600   
 #     -1963//204120     
 # 226287557//37623398400

macro convfrac(T,n,d)
    return quote
        convert($(esc(T)),$(esc(n))) / convert($(esc(T)),$(esc(d)))
    end
end
      
# Coefficients calculated by hand
# const before constants slows exection
function wser(p,ps)
    T = typeof(p)
    μ₂ = one(T)/3   # minus mu2 !
    μ₃ = @convfrac(T,11,72) # must make var to get compiler optimization (v0.3)
    μ₄ = @convfrac(T,-43,540)
    μ₅ = @convfrac(T,769,17280)
    μ₆ = @convfrac(T,680863,43545600)
    p - ps * (μ₂ - p * (μ₃ + p * (μ₄ + p * (μ₅ + p * μ₆))))
end

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

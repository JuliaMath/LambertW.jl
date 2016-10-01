using DeepConvert
# lambertw loses precision near -1/e
# These are diagnostics.
# One thing that at least helps convergence near the pole
# is to add a very small imaginary part.  1e-100 * im or so.

function converr(z,k)
    w = lambertw(z,k)
   convert(Float64,abs(w*exp(w)-z))
end

converr(z) = converr(z,0)

mbig{T<:Real}(x::T) = BigFloat(x)
mbig(z::Complex) = complex(BigFloat(real(z)),BigFloat(imag(z)))

# compare lambertw to bigfloat version
function bigcomp(z,k)
    zr = eval(parse(z))
    zb = deepbigfloat(z)
    wb = lambertw(zb,k)
    w = lambertw(zr,k)
    -log10(convert(Float64,abs(wb-w)))
end


bigcomp(z) = bigcomp(z,0)

# compare lambertw to series expansion about pole
function sercomp(x::String)
    xr = eval(parse(x))
    xb = deepbigfloat(x) - 1/big(e)
    w = lambertw(xb,-1)
    ws = lambertwm1(xr)
    -log10(convert(Float64,abs(ws-w)))
end

# compare sercomp("1/10^7") to bigcomp("-1/e + 1/10^7",-1)
# They are roughly the same. So this is the crossover point:
# x = 1/10^7
# for smaller x, the series lambertwm1 is better for k=-1.

using LambertW

mu1 = 1
mu2 = -1//3
mu3 = 1//8 * (-2*mu2 + mu1^2 -4 * mu2^2)
mu4 = 1//10 * (-2*mu3 + 3 * mu1 * mu2 -10 * mu2 * mu3)
mu5 = 1//12 * (-2*mu4 + 2 * (2*mu1*mu3+mu2^2) - 6 * (mu3^2 + 2*mu2*mu4))

function muterm(n,a)
    T = eltype(a)
    sum1 = zero(T)
    for j in 1:n-1
        sum1 += a[j] * a[n-j]
#        println("sum1 $sum1")
    end
    sum1 *= n/convert(T,2)
#    println("  * sum1 $sum1")
    sum2 = zero(T)
    for j in 2:n
        sum2 += a[j] * a[n+2-j]
#        println("sum2 $sum2")
    end
    sum2 *= -(n+2)
#    println("  * sum2 $sum2")
    a[n+1] = (-2*a[n] + sum1 + sum2)/(2*(n+2))
end

function genmu(T::DataType, n::Int)
    a = Array(T,n)
    a[1] =  one(T)
    a[2] = -one(T) / convert(T,3)
#    a[3] = convert(T,11) / convert(T,72)
    for k in 2:n-1
        muterm(k,a)
    end
    a
end

mu64a = genmu(Float64,6)
mu64b = genmu(Rational,6)
mu64c = genmu(BigFloat,100)

# lambertwbf with many terms
function newlam(x)
    a = mu64a
    s = zero(x)
    y = sqrt(2*e*x)
#    y = 1.0
    for i in 1:length(a)
#        println("mu $(a[i])")
        s += y^i * a[i]
    end
    s
end

function hornerlam(x)
    a = mu64a
    y = sqrt(2*e*x)
    s = zero(x)
    for i in (length(a)):-1:1
        s = a[i] + s * y
    end
    return s * y
end


# Test shows backward is more precise than forwards.
function rnewlam(x)
    s = zero(x)
    y = sqrt(2*e*x)
    for i in length(mu64a):-1:1
        s += y^i * mu64a[i]
    end
    s
end

function bnewlam(x)
    s = zero(x)
    y = sqrt(2*e*x)
    for i in 1:length(mu64c)
        s += y^i * mu64c[i]
    end
    s
end


function lamnew(y)
    x = sqrt(2*e*y)
    mu1*x + mu2*x^2 + mu3*x^3 + mu4*x^4 +mu5*x^5
end

lamwex(x) = 1+lambertw(-1/big(e)+BigFloat(1)/x)

lamwex64(x) = 1+lambertw(-1/e + 1/x)

diffx(x) = float64((w=lamwex(x); (lambertwbp(float64(1/x)) - w)/w))
diffxn(x) = float64((w=lamwex(x); (lamnew(float64(1/x)) - w)/w))

diffrat(x) = abs(diffxn(x)/diffx(x))

diffx64(x) = (w=lamwex(x); float64((w-lamwex64(x))/w))

diffrat64n(x) = abs(diffxn(x)/diffx64(x))
diffrat64(x) = abs(diffx(x)/diffx64(x))

# 1. Compute the "exact" value of lambertw/
# 2. Compute series values with coeffs from paper, and our coeffs.
# 3. Compute error in each series value.
# 4. Print ratio of errors.
# Small number means new coeffs are doing better.
# column 2 is x = z + 1/e.
# column 3 is the value for the variable in series expansion: sqrt(2*e*x)
# column 4 is ratio of error in new series to errr in old series
# column 5 is ratio of error in old series to exact val Float64 error
# column 6 is ratio of error in new series to exact val Float64 error
# Crossover point for using series with Float64 is where column 5 is about unity
# This is about x = 1e-5
# Results show increasing improvement until higher order terms are not
# significant.
# We use inverses here and there to avoid losing BigFloat precision
function diffratpr()
    println("y        x        p         en/eo     eo/eex    en/eex")
    iex = 0
    pref = 2
    for i in 1:30
        if iseven(i)
            pref = 7
        else
            iex += 1
            pref = 1
        end
        y = pref*10^iex
        x = 1/y
        @printf( "%.1e  %.1e  %.2e  %.2e  %.2e  %.2e\n" ,
                float(y), x , sqrt(2*e*x),  diffrat(y), diffrat64(y),diffrat64n(y))
    end
end



### Following is only playing with rewriting expressions.

iscall(ex::Expr) = ex.head == :call
isop(ex::Expr, s::Symbol) = iscall(ex) && ex.args[1] == s
function is_x_times_fofx(ex::Expr, f::Symbol)
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

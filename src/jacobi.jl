export jacobisymbol

macro flip!(n) :($(esc(n)) = -$(esc(n))) end
function jacobisymbol(a::Union(Signed,Unsigned),n::Union(Signed,Unsigned))
    if n <= 0 || iseven(n) throw(DomainError()) end    
    j = 1
    if n < 0
        n % 4 == 3  ? @flip!(j)  : nothing
        @flip!(n)
    end
    while a != 0
        while iseven(a)
            n % 8 == 3 || n % 8 == 5 ? @flip!(j) : nothing
            a >>= 1            
        end
        (a,n) = (n,a)
        a % 4 == 3 && n % 4 == 3 ? @flip!(j) : nothing
        a %= n
    end    
    return n == 1 ? j : 0
end

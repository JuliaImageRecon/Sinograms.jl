function outer_sum(xx, yy)
#=
function ss = outer_sum(xx,yy)

 compute an "outer sum" x + y'
 that is analogous to the "outer product" x * y'

in
	xx	[nx 1]
	yy	[1 ny]
		more generally: xx [(dim)] + yy [L,1] -> xx [(dim) LL]
out
	ss [nx ny]	ss(i,j) = xx(i) + yy(j)

Translated from outer_sum.m in MIRT
Copyright 2022, Jeff Fessler and Jason Hu, University of Michigan
=#
    xx = xx[:]
    yy = yy[:]
    nx = length(xx)
    ny = length(yy)
    ans = zeros(nx,ny)
    for i = 1:nx
        for j = 1:ny
            ans[i,j] = xx[i]+yy[j]
        end
    end
    return ans
end

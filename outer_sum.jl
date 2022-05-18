function outer_sum(xx, yy)
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

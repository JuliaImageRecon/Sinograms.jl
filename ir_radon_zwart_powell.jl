## Matlab author: Seongjin Yoon, Univ. of Michigan, Jan 2015
## Translated by: Harshit Nanda

### Analytic 2D Radon transform of Zwart-Powell box spline.

#=
Reference
[1] A. Entezari, M. Nilchian, and M. Unser, "A box spline calculus
for the discretization of computed tomography reconstruction problems,"
IEEE Trans. Med. Imaging, vol. 31, no. 8, pp. 1532–1541, 2012.
2015-08-10 Jeff Fessler, added self test and parallelized

## Mathematcatical representation

# theta: ray angle in radian
# rr: distance between the point and the ray (normalized by the pixel size)
# output: radon transform of the Zwart-Powell box spline element

=#


function ir_radon_zwart_powell(theta, rr)
  dim = size(theta)
  theta = vec(theta)    
  zeta = zeros(length(theta), 4) 
  zeta[:,1] = cos.(theta)
  zeta[:,2] = sin.(theta)
  zeta[:,3] .= zeta[:,1] .+ zeta[:,2]
  zeta[:,4] .= zeta[:,1] .- zeta[:,2] 
  cond = abs.(zeta) .>= 1.1920929e-07  #eps('single') matlab different than eps(1.0) in julia
  N = sum(cond, dims = 2)
  N = vec(N)      
  output = BoxSp4(rr[:], zeta, cond, N) ./ factorial.(N.-1)
  output = reshape(output, dim)
  return output
end 
    
BoxSp0(y, N::Int) = y > 0 ? y^(N-1) : zero(y)

#=
function output = BoxSp0(y, N)
  output = heaviside(y) .* y.^(N-1);
  if any(size(y) ~= size(N)), 
    keyboard, 
  end
  output = (y >= 0) .* y.^(N-1);
end 
=#
    

function BoxSp1(y, zeta, cond, N)
  good = cond[:,1]
  output = (BoxSp0.(y .+ (0.5 .* zeta[:,1]), N) .- BoxSp0.(y .- (0.5 .* zeta[:,1]), N)) ./ zeta[:,1]
  output[.~good] .= BoxSp0.(y[.~good], N[.~good])
  return output
end

#= 
function output = BoxSp1(y, zeta, cond, N)
  good = cond(:,1);
  output = (BoxSp0(y+0.5*zeta(:,1), N) - BoxSp0(y-0.5*zeta(:,1), N)) ./ zeta(:,1);
  output(~good) = BoxSp0(y(~good), N(~good));
end 
=#
    
function BoxSp2(y, zeta, cond, N)
  good = cond[:,2]
  output = (BoxSp1(y .+ (0.5.*zeta[:,2]), zeta, cond, N) .- BoxSp1(y .- (0.5 .* zeta[:,2]), zeta, cond, N)) ./ zeta[:,2]
  output[.~good] .= BoxSp1(y[.~good], zeta[.~good,:], cond[.~good,:], N[.~good])
  return output 
end

#= 
function output = BoxSp2(y, zeta, cond, N)
  good = cond(:,2);
  output = (BoxSp1(y+0.5*zeta(:,2), zeta, cond, N) - BoxSp1(y-0.5*zeta(:,2), zeta, cond, N)) ./ zeta(:,2);
  output(~good) = BoxSp1(y(~good), zeta(~good,:), cond(~good,:), N(~good));
end 
=#
    
function BoxSp3(y, zeta, cond, N)
  good = cond[:,3]
  output = (BoxSp2(y .+ 0.5 .* zeta[:,3], zeta, cond, N) .- BoxSp2(y .- 0.5 .* zeta[:,3], zeta, cond, N)) ./ zeta[:,3]
  output[.~good] .= BoxSp2(y[.~good], zeta[.~good,:], cond[.~good,:], N[.~good])
  return output
end

#=
function output = BoxSp3(y, zeta, cond, N)
  good = cond(:,3);
  output = (BoxSp2(y+0.5*zeta(:,3), zeta, cond, N) - BoxSp2(y-0.5*zeta(:,3), zeta, cond, N)) ./ zeta(:,3);
  t = zeta(~good,:);
  output(~good) = BoxSp2(y(~good), t, cond(~good,:), N(~good));
end 
=#
    
function BoxSp4(y, zeta, cond, N)
  good = cond[:,4]
  output = (BoxSp3(y .+ 0.5 .* zeta[:,4], zeta, cond, N) .- BoxSp3(y .- 0.5 .* zeta[:,4], zeta, cond, N)) ./ zeta[:,4] 
  output[.~good] .= BoxSp3(y[.~good], zeta[.~good,:], cond[.~good,:], N[.~good])
  return output
end

#=
function output = BoxSp4(y, zeta, cond, N)
  good = cond(:,4);
  output = (BoxSp3(y+0.5*zeta(:,4), zeta, cond, N) - BoxSp3(y-0.5*zeta(:,4), zeta, cond, N)) ./ zeta(:,4);
  output(~good) = BoxSp3(y(~good), zeta(~good,:), cond(~good,:), N(~good));
end =#
    
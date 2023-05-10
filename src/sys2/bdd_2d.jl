## Author: Sonia Minseo Kim
## Reference: Rodrigo de Barros Vimieiro

## 2-D Branchless Distance Driven Code
using LazyGrids: ndgrid 
using Interpolations: linear_interpolation, Flat 

## Projection
"""
    projection(phantom::AbstractMatrix{<:T}, geo)
Generates a sinogram using the forward projection algorithm,
for a phantom image and a tuple of geometry definitions 
"""
function projection(phantom::AbstractMatrix{<:T}, geo ; draw::Bool = false) where {T <: Number}

    DSD = geo.DSD
    DS0 = geo.DS0      
    pSize = geo.pSize    
    dSize = geo.dSize    
    nPix = geo.nPix      
    nDet = geo.nDet     
    angle = geo.angle # vector of projection view angles in radians

    # Detector boundaries
    detX = (-(nDet/2):(nDet/2)) .* dSize
    detY = (-(DSD-DS0)-(dSize/2)) .* ones(nDet+1)

    # Pixel boundaries
    (pixelX,pixelY) = ndgrid((-(nPix/2):(nPix/2)), (-(nPix/2):(nPix/2)))
    pixelX = pixelX' * pSize
    pixelY = reverse(pixelY' .- 1/2, dims = 1)
    pixelY = pixelY * pSize

    sinogram = zeros(length(angle),nDet)

    # For each projection
    for proj in 1:length(angle)
        
        beta = angle[proj] # angle from x-ray beam to y-axis
        
        # Tube rotation
        rtubeX = -DS0*sin(beta)
        rtubeY = DS0*cos(beta)

        # Detector rotation
        rdetX = detX.*cos(beta) - detY.*sin(beta)
        rdetY = detX.*sin(beta) + detY.*cos(beta)

        # Define angle case & which axis it project boundaries
        if (((beta >= π/4) && (beta <= 3*π/4)) || ((beta >= 5*π/4) && (beta <= 7*π/4)))
            axisCase = false # map on y axis
        else
            axisCase = true # map on x axis
        end
        
        # Mapping boundaries onto a common axis
        if (axisCase)
            detm = map2x.(rtubeX,rtubeY,rdetX,rdetY)
            pixm = map2x.(rtubeX,rtubeY,pixelX,pixelY)
            img = phantom
        else
            detm = map2y.(rtubeX,rtubeY,rdetX,rdetY)
            pixm = reverse(map2y.(rtubeX,rtubeY,pixelX,pixelY)', dims = 2)
            img = reverse(phantom', dims = 2)
        end

        # scaling factor calculation
        detSize = diff(detm, dims = 1)
        L = zeros(1,nDet) 

        if (axisCase) 
            for n = 1:nDet 
                det_mid = (detm[n] + detm[n+1])/2
                theta = atan(abs(rtubeX-det_mid)/abs(rtubeY)) # theta is the angle btw the ray to det_mid and y-axis
                L[n] = abs(pSize./(detSize[n].*cos(theta)))
            end
        else
            for n = 1:nDet
                det_mid = (detm[n] + detm[n+1])/2
                theta = acot(abs(rtubeY-det_mid)/abs(rtubeX))
                L[n] = abs(pSize./(detSize[n].*sin(theta)))
            end 
        end     
        
        sinoTmp = zeros(1, nDet) # one "row" of sinogram

        # For each image row
        for row in 1:nPix
            
            rowm = pixm[row,:] # mapped row from pixel mapped         
            pixSize = diff(rowm, dims = 1)

            Ppj = integrate1D(img[row,:], pixSize)
            
            interp_func = linear_interpolation(rowm, Ppj, extrapolation_bc=Flat())
            Pdk = interp_func.(detm)

            sinoTmp = sinoTmp + abs.(diff(Pdk, dims = 1)');              
            
        end # Row loop 
        sinogram[proj,:] = sinoTmp .* L

    end # Projection loop

    return sinogram

end # endfunc


## Backprojection 
"""
    backprojection(sinogram::AbstractMatrix{<:T}, geo)
Generates a reconstructed image using the back projection algorithm,
for a sinogram and a tuple of geometry definitions 
"""
function backprojection(sinogram::AbstractMatrix{<:T},geo ; draw::Bool = false) where {T <: Number}

    DSD = geo.DSD      
    DS0 = geo.DS0      
    pSize = geo.pSize    
    dSize = geo.dSize    
    nPix = geo.nPix       
    nDet = geo.nDet     
    angle = geo.angle 

    # Detector boundaries
    detX = (-(nDet/2):(nDet/2)) .* dSize
    detY = (-(DSD-DS0)-(dSize/2)) .* ones(nDet+1)

    # Pixel boundaries
    (pixelX,pixelY) = ndgrid(-(nPix/2):(nPix/2), -(nPix/2):(nPix/2))
    pixelX = pixelX' * pSize
    pixelY = reverse(pixelY' .- 1/2, dims = 1)
    pixelY = pixelY * pSize

    reconImg = zeros(nPix,nPix)
    reconImgTmp = reconImg

    # For each projection
    for proj in 1:length(angle)
        
        beta = angle[proj]

        # Tube rotation
        rtubeX = -DS0*sin(beta)
        rtubeY = DS0*cos(beta)

        # Detector rotation
        rdetX = detX.*cos(beta) - detY.*sin(beta)
        rdetY = detX.*sin(beta) + detY.*cos(beta)

        # Define angle case & which axis it project boundaries
        if (((beta >= π/4) && (beta <= 3*π/4)) || ((beta >= 5*π/4) && (beta <= 7*π/4)))
            axisCase = false # map on y axis
        else
            axisCase = true # map on x axis
        end
        
        # Mapping boundaries onto a common axis
        if (axisCase)
            detm = map2x.(rtubeX,rtubeY,rdetX,rdetY)
            pixm = map2x.(rtubeX,rtubeY,pixelX,pixelY)
        else
            detm = map2y.(rtubeX,rtubeY,rdetX,rdetY)
            pixm = reverse(map2y.(rtubeX,rtubeY,pixelX,pixelY)', dims = 2)
        end
        
        # scaling factor calculation
        L = zeros(1,nPix) 
        if (axisCase) 
            for n = 1:nPix 
                pixSize = pixm[n,2]-pixm[n,1] 
                pix_mid = (pixm[n] + pixm[n+1])/2
                theta = atan(abs(rtubeX-pix_mid)/abs(rtubeY)) # theta is the angle btw the ray to det_mid and y-axis
                L[n] = abs(pSize./(pixSize.*cos(theta)))
            end
        else
            for n = 1:nPix
                pixSize = pixm[n,2]-pixm[n,1] 
                pix_mid = (pixm[n] + pixm[n+1])/2
                theta = acot(abs(rtubeY-pix_mid)/abs(rtubeX))
                L[n] = abs(pSize./(pixSize.*sin(theta)))
            end  
        end     

        # For each row
        for row in 1:nPix
            
            rowm = pixm[row,:] # Get first mapped row from image --> changed to column vec
            detSize = diff(detm, dims = 1)      
            
            #detm = vec(detm)
            Ppj = integrate1D(sinogram[proj,:],detSize)

            if detm[1] > detm[2] #if descending
                interp_func = linear_interpolation(reverse(detm), reverse(vec(Ppj)), extrapolation_bc=Flat())
            else
                interp_func = linear_interpolation(detm, vec(Ppj), extrapolation_bc=Flat())
            end
            Pdk = interp_func.(vec(rowm))
            
            reconImgTmp[row,:] = abs.(diff(Pdk, dims = 1)') .* L
            
        end # Row loop 
        
        if (axisCase)
            reconImg = reconImg + reconImgTmp 
        else
            reconImg = reconImg + reverse(reconImgTmp, dims = 2)'
        end 
           
    end # Projection loop
    
    reconImg = reconImg / length(angle) 
    return reconImg

end # endfunc


## Digital Integration function
"""
    integrate1D(p_v::Vector, pixelSize)
Calculates the integral image set,
for a given column vector and a pixel size 
"""
function integrate1D(p_v::Vector,pixelSize::Vector)
    n_pixel = length(p_v)
    
    P_x = 0
    Ppj = zeros(n_pixel+1)
    
    for pj in 1:n_pixel
       P_x += p_v[pj] * pixelSize[pj]
       
       Ppj[pj+1] = P_x
    end
    return Ppj
end

## Map Y
"""
    map2y(x1,y1,x2,y2)
Maps detector or pixel boundaries onto y-axis,
for tube and detector rotation angles and detector/pixel boundaries
"""
function map2y(x1,y1,x2,y2)
    y = y1-x1*(y1-y2)/(x1-x2)
    return y
end

## Map X
"""
    map2x(x1,y1,x2,y2)
Maps detector or pixel boundaries onto x-axis,
for tube and detector rotation angles and detector/pixel boundaries
"""
function map2x(x1,y1,x2,y2)
    x = x1-y1*(x1-x2)/(y1-y2)
    return x
end
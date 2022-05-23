## Matlab Author: Rodrigo de Barros Vimieiro
## Translated by Sonia Minseo Kim in 2022

## 2-D Branchless Distance Driven Code

export bdd_2d

# using Plots: plot
using ImagePhantoms: shepp_logan, SheppLoganToft
using LazyGrids: ndgrid 
using Dierckx: Spline1D

## Projection
"""
    projection(phantom::AbstractMatrix{<:T}, geo)
Generates a sinogram using the forward projection algorithm,
for a phantom image and a tuple of geometry definitions 
"""
function projection(phantom::AbstractMatrix{<:T}, geo ; draw::Bool = false) where {T <: Number}

    DSD = geo.DSD
    DSO = geo.DSO      
    pSize = geo.pSize    
    dSize = geo.dSize    
    nPix = geo.nPix      
    nDet = geo.nDet     
    theta = geo.theta # vector of projection view angles in radians

    # Detector boundaries
    detX = (-(nDet/2):(nDet/2)) .* dSize
    detY = (-(DSD-DSO)-(dSize/2)) .* ones(nDet+1)

    # Pixel boundaries
    (pixelX,pixelY) = ndgrid((-(nPix/2):(nPix/2)), (-(nPix/2):(nPix/2)))
    pixelX = pixelX' .* pSize
    pixelY = reverse(pixelY' .- pSize/2, dims = 1)

    # Tube
    tubeX = 0     
    tubeY = DSO    

    # Iso-center
    isoX = geo.isoX
    isoY = geo.isoY

    sinogram = zeros(length(theta),nDet)

    # For each projection
    for proj in 1:length(theta)
        
        angle = theta[proj]
        
        # Tube rotation
        rtubeX = ( (tubeX - isoX )*cos(angle) - (tubeY - isoY )*sin(angle) ) + isoX
        rtubeY = ( (tubeX - isoX )*sin(angle) + (tubeY - isoY )*cos(angle) ) + isoY

        # Detector rotation
        rdetX = ( (detX .- isoX ).*cos(angle) - (detY .- isoY ).*sin(angle) ) .+ isoX
        rdetY = ( (detX .- isoX ).*sin(angle) + (detY .- isoY ).*cos(angle) ) .+ isoY
        
        if (draw == true)
            # drawGeo(rtubeX,rtubeY,rdetX,rdetY)   
        end

        # Define angle case & which axis it project boundaries
        # Case 1 
        if (((angle>=0)&&(angle<=pi/4))||(angle>=7*pi/4))     
            axisXCase = true   # Map on X axis()
            angleCase = 1
            c1=0
            c2=1
        # Case 2
        elseif ((angle>pi/4) && (angle<3*pi/4))
            axisXCase = false   # Map on Y axis()
            angleCase = 2
            c1=0
            c2=1
        # Case 3
        elseif ((angle>=3*pi/4)&&(angle<=5*pi/4))
            axisXCase = true   # Map on X axis()
            angleCase = 3
            c1=0
            c2=-1
        else # Case 4
            axisXCase = false   # Map on Y axis()
            angleCase = 4
            c1=0
            c2=-1
        end
        
        # Mapping boundaries into a common axis()
        if (axisXCase)
            detm = mapp2x.(rtubeX,rtubeY,rdetX,rdetY)
            pixm = mapp2x.(rtubeX,rtubeY,pixelX,pixelY)
            img = phantom
        else
            detm = mapp2y.(rtubeX,rtubeY,rdetX,rdetY)
            pixm = reverse(mapp2y.(rtubeX,rtubeY,pixelX,pixelY)', dims = 2)
            img = reverse(phantom', dims = 2)
        end
                
        center_det = floor(Int, (nDet+1)/2+1)

        if (axisXCase) 
            # X-Ray pixel intersection calculation
            L = abs(pSize/cos(angle)) # This account for parallel-beam
            # Correction for fan-beam
            L1 = @. sqrt((rtubeX-detm[center_det])^2+(rtubeY-0)^2)/sqrt((rtubeX-detm[1:nDet]')^2+(rtubeY-0)^2)
        else
            # X-Ray pixel intersection calculation
            L = abs(pSize/sin(angle)) # This account for parallel-beam
            # Correction for fan-beam
            L1 = @. sqrt((rtubeX-0)^2+(rtubeY-detm[center_det])^2)/sqrt((rtubeX-0)^2+(rtubeY-detm[1:nDet]')^2)       
        end     
        L = L ./ L1
        
        pixIstart = 1
        pixIinc = 1
        if ((angleCase == 1)||(angleCase == 2))
            detIstart = 1
            detIinc = 1
        else
            detIstart = nDet+1
            detIinc = -1
        end

        
        deltaDetm = detm[detIstart+detIinc] - detm[detIstart] # Mapped detector length()
        deltaPixm = pixm[1,2] - pixm[1,1]   # Mapped pixel length()
        
        sinoTmp = zeros(1, nDet) # one "row" of sinogram

        # For each image row
        for row in 1:nPix
            
            rowm = pixm[row,:] # Get first mapped row from image

            detInd = detIstart
            pixInd = pixIstart            
            
            # Find first detector overlap mapped with pixel mapped [Case 1]
            if (detm[detInd]-rowm[pixIstart]<=deltaDetm)
                while detm[detInd]-rowm[pixIstart]<=deltaDetm            
                    detInd = detInd + detIinc            
                end
            else
            # Find first pixel overlap mapped with detector mapped [Case 2]           
                if (detm[detIstart]-rowm[pixInd]>deltaPixm)
                    while detm[detIstart]-rowm[pixInd]>deltaPixm            
                        pixInd = pixInd + pixIinc            
                    end
                end
            end
            
            Ppj = integrate1D(img[row,:],deltaPixm)

            if rowm[1] > rowm[2] #if descending
                Pdk = Spline1D(reverse(rowm), reverse(Ppj); k=1, bc="zero")(detm)
                Pdk = reverse(Pdk)
            else 
                Pdk = Spline1D(rowm, Ppj; k=1, bc="zero")(detm)
            end
            
            while detm[detInd+c2]<rowm[end]
                sinoTmp[detInd] = sinoTmp[detInd]+(Pdk[detInd+c2]-Pdk[detInd])/deltaDetm
                detInd = detInd + detIinc
            end
            sinoTmp[detInd] = sinoTmp[detInd]+(Ppj[end]-Pdk[detInd])/deltaDetm                   
            
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
    DSO = geo.DSO      
    pSize = geo.pSize    
    dSize = geo.dSize    
    nPix = geo.nPix       
    nDet = geo.nDet     
    theta = geo.theta 

    # Detector boundaries
    detX = (-(nDet/2):(nDet/2)) .* dSize
    detY = (-(DSD-DSO)-(dSize/2)) .* ones(nDet+1)

    # Pixel boundaries
    (pixelX,pixelY) = ndgrid(-(nPix/2):(nPix/2), -(nPix/2):(nPix/2))
    pixelX = pixelX' .* pSize
    pixelY = reverse(pixelY' .- pSize/2, dims = 1)

    # Tube
    tubeX = 0     
    tubeY = DSO    

    # Iso-center
    isoX = geo.isoX
    isoY = geo.isoY

    reconImg = zeros(nPix,nPix)
    reconImgTmp = reconImg

    # For each projection
    for proj in 1:length(theta)
        
        angle = theta[proj]

        # Tubre rotation
        rtubeX = ( (tubeX - isoX )*cos(angle) - (tubeY - isoY )*sin(angle) ) + isoX
        rtubeY = ( (tubeX - isoX )*sin(angle) + (tubeY - isoY )*cos(angle) ) + isoX

        # Detector rotation
        rdetX = ( (detX .- isoX )*cos(angle) - (detY .- isoY )*sin(angle) ) .+ isoX
        rdetY = ( (detX .- isoX )*sin(angle) + (detY .- isoY )*cos(angle) ) .+ isoX
        
        if (draw == true)
            # drawGeo(rtubeX,rtubeY,rdetX,rdetY)
        end    
        
        # Define angle case & which axis it project boundaries
        # Case 1 
        if (((angle>=0)&&(angle<=pi/4))||(angle>=7*pi/4))       
            axisXCase = true   # Map on X axis()
            angleCase = 1
            c1=1
            c2=0
        # Case 2
        elseif ((angle>pi/4)&&(angle<3*pi/4))
            axisXCase = false   # Map on Y axis()
            angleCase = 2
            c1=1
            c2=0
        # Case 3
        elseif (((angle>=3*pi/4)&&(angle<=5*pi/4)))
            axisXCase = true   # Map on X axis()
            angleCase = 3
            c1=0
            c2=1
        else # Case 4
            axisXCase = false   # Map on Y axis()
            angleCase = 4
            c1=0
            c2=1
        end
        
        # Mapping boundaries into a common axis()
        if (axisXCase)
            detm = mapp2x.(rtubeX,rtubeY,rdetX,rdetY)
            pixm = mapp2x.(rtubeX,rtubeY,pixelX,pixelY)
        else
            detm = mapp2y.(rtubeX,rtubeY,rdetX,rdetY)
            pixm = reverse(mapp2y.(rtubeX,rtubeY,pixelX,pixelY)', dims = 2)
        end
        
                
        center_det = floor(Int, (nPix+1)/2+1)
        
        if (axisXCase) 
            # X-Ray pixel intersection calculation
            L = abs(pSize/cos(angle)) # This account for parallel-beam
            # Correction for fan-beam
            L1 = @. sqrt((rtubeX-pixm[center_det])^2+(rtubeY-0)^2)/sqrt((rtubeX-pixm[1:nPix]')^2+(rtubeY-0)^2)
        else
            # X-Ray pixel intersection calculation
            L = abs(pSize/sin(angle)) # This account for parallel-beam
            # Correction for fan-beam
            L1 = @. sqrt((rtubeX-0)^2+(rtubeY-pixm[center_det])^2)/sqrt((rtubeX-0)^2+(rtubeY-pixm[1:nPix]')^2)          
        end     
        L = L ./ L1
        
        pixIstart = 1
        pixIinc = 1
        if ((angleCase == 1)||(angleCase == 2))
            detIstart = 1
            detIinc = 1
        else
            detIstart = nDet+1
            detIinc = -1
        end

        
        deltaDetm = detm[detIstart+detIinc]- detm[detIstart] # Mapped detector length()
        deltaPixm = pixm[1,2]- pixm[1,1]   # Mapped pixel length()    
        

        # For each row
        for row in 1:nPix
            
            rowm = pixm[row,:] # Get first mapped row from image --> changed to column vec
            reconTmp = zeros(1,nPix)

            detInd = detIstart
            pixInd = pixIstart            
            
            # Find first detector overlap mapped with pixel mapped [Case 1]
            if (detm[detInd]-rowm[pixIstart]<=deltaDetm)
                while ((detm[detInd]-rowm[pixIstart]<=deltaDetm))            
                    detInd = detInd + detIinc            
                end
            else
            # Find first pixel overlap mapped with detector mapped [Case 2]           
                if (detm[detIstart]-rowm[pixInd]>deltaPixm)
                    while (detm[detIstart]-rowm[pixInd]>deltaPixm)            
                        pixInd = pixInd + pixIinc           
                    end
                end
            end
            
            detm = vec(detm)
            Ppj = integrate1D(sinogram[proj,:],deltaDetm)

            if detm[1] > detm[2] #if descending
                Pdk = Spline1D(reverse(detm), reverse(vec(Ppj)); k=1, bc="zero")(vec(rowm))
            else 
                Pdk = Spline1D(detm, vec(Ppj); k=1, bc="zero")(vec(rowm))
            end
            
            while (pixInd < nPix+1)
                reconTmp[pixInd] = reconTmp[pixInd]+(Pdk[pixInd+c1]-Pdk[pixInd+c2])/deltaDetm
                pixInd = pixInd + pixIinc 
            end
            
            reconImgTmp[row,:] = reconTmp .* L
            
        end # Row loop 
        
        if ((angleCase == 1)||(angleCase == 3))
            reconImg = reconImg + reconImgTmp 
        else
            reconImg = reconImg + reverse(reconImgTmp, dims = 2)'
        end 
           
    end # Projection loop
    
    reconImg = reconImg / length(theta) 
    return reconImg

end # endfunc



## Integrate function
"""
    integrate1D(p_v::Vector, pixelSize)
Calculates the integral image set,
for a given column vector and a pixel size 
"""
function integrate1D(p_v::Vector,pixelSize)

    n_pixel = length(p_v)
    
    P_x = 0
    Ppj = zeros(n_pixel+1)
    
    for pj in 1:n_pixel

       P_x += p_v[pj] * pixelSize
       
       Ppj[pj+1] = P_x
       
    end

    return Ppj

end

## Map Y
"""
    mapp2y(x1,y1,x2,y2)
Maps detector or pixel boundaries onto y-axis,
for tube and detector rotation angles and detector/pixel boundaries
"""
function mapp2y(x1,y1,x2,y2)
    y = y1-x1*(y1-y2)/(x1-x2)
    return y
end

## Map X
"""
    mapp2x(x1,y1,x2,y2)
Maps detector or pixel boundaries onto x-axis,
for tube and detector rotation angles and detector/pixel boundaries
"""
function mapp2x(x1,y1,x2,y2)
    x = -(y1)*(x1-x2)/(y1-y2)+x1
    return x
end

## Draw 2D
# Function to draw geometry
#=
function drawGeo(tubeX,tubeY,detX,detY)
    plot(tubeX,tubeY,markershape = :asterisk)
    plot!([detX[1],detX[end]],[detY[1],detY[end]])
end
=#

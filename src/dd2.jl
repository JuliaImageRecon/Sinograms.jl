## Matlab Author: Rodrigo de Barros Vimieiro
## Translated by Sonia Minseo Kim

## 2-D Branchless Distance Driven Code
using Plots: plot
using ImagePhantoms: shepp_logan, SheppLoganToft
using ImageGeoms: ImageGeom, axesf
using MIRTjim: jim, prompt #specify the functions you are using in the packages
using LazyGrids: ndgrid 
using Dierckx: Spline1D

# Julia is compiling program so you need to restart the program in order to clear the variables 

## Geometry Definitions
# Geometry
deg = 1
geo = (DSD = 100000, DSO = 99700, pSize = 1, dSize = 0.5, nPix = 256, nDet = 1024
        , theta = deg2rad.(0:deg:360-deg), isoX = 0, isoY = 0) #Define tuple

# Make the Phantom Image
phantomImg = shepp_logan(geo.nPix, SheppLoganToft())

## Projection Branchless
function projectionBranchless(phantom,geo ; draw::Bool = false)

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
pixelX = pixelX .* pSize
pixelY = reverse(pixelY .- pSize/2, dims = 1)

# Tube
tubeX = 0     
tubeY = DSO    

# Iso-center
isoX = geo.isoX
isoY = geo.isoY

sinogram = zeros(size(theta,1),nDet)

# For each projection
for proj=1:size(theta,1)
    
    angle = theta[proj]
       
    # Tubre rotation
    rtubeX = ( (tubeX - isoX )*cos(angle) - (tubeY - isoY )*sin(angle) ) + isoX
    rtubeY = ( (tubeX - isoX )*sin(angle) + (tubeY - isoY )*cos(angle) ) + isoY

    # Detector rotation
    rdetX = ( (detX .- isoX ).*cos(angle) - (detY .- isoY ).*sin(angle) ) .+ isoX
    rdetY = ( (detX .- isoX ).*sin(angle) + (detY .- isoY ).*cos(angle) ) .+ isoY
    
    if (draw == true)
        drawGeo(rtubeX,rtubeY,rdetX,rdetY)   
    end

    # Define angle case & which axis it it project boundaries
    # Case 1 
    if (((angle>=0)&&(angle<=pi/4))||(angle>=7*pi/4))     
        axisXCase = true   # Map on X axis()
        angleCase = 1
        c1=0
        c2=1
    else
        # Case 2
        if ((angle>pi/4) && (angle<3*pi/4))
            axisXCase = false   # Map on Y axis()
            angleCase = 2
            c1=0
            c2=1
        else
            # Case 3
            if ((angle>=3*pi/4)&&(angle<=5*pi/4))
                axisXCase = true   # Map on X axis()
                angleCase = 3
                c1=0
                c2=-1
            else
            # Case 4
                axisXCase = false   # Map on Y axis()
                angleCase = 4
                c1=0
                c2=-1
            end
        end
    end
    
    # Mapping boundaries into a commum axis()
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
    L1 = zeros(1, nDet) 

    if (axisXCase) 
        # X-Ray pixel intersection calculation
        L = abs(pSize/cos(angle)) # This account for parallel-beam
        # Correction for fan-beam
        for n=1:nDet
            L1[n]=sqrt((rtubeX-detm[center_det])^2+(rtubeY-0)^2)/sqrt((rtubeX-detm[n])^2+(rtubeY-0)^2)
        end
    else
        # X-Ray pixel intersection calculation
        L = abs(pSize/sin(angle)) # This account for parallel-beam
        # Correction for fan-beam
        for n=1:nDet
            L1[n]=sqrt((rtubeX-0)^2+(rtubeY-detm[center_det])^2)/sqrt((rtubeX-0)^2+(rtubeY-detm[n])^2)
        end           
    end     
    L = L./L1
    
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
    
    sinoTmp = zeros(1,nDet)
    
    # For each row
    for row in 1:nPix
        
        rowm = pixm[row,:] # Get first mapped row from image.
        
        detInd = detIstart
        pixInd = pixIstart            
        
        # Find first detector overlap maped with pixel maped [Case 1]
        if (detm[detInd]-rowm[pixIstart]<=deltaDetm)
            while detm[detInd]-rowm[pixIstart]<=deltaDetm            
                detInd = detInd + detIinc            
            end
        else
        # Find first pixel overlap maped with detector maped [Case 2]           
            if (detm[detIstart]-rowm[pixInd]>deltaPixm)
                while detm[detIstart]-rowm[pixInd]>deltaPixm            
                    pixInd = pixInd + pixIinc            
                end
            end
        end

        Ppj = integrate1D(img[row,:],deltaPixm)
        Pdk = Spline1D(sort!(rowm), vec(Ppj); k=1)(detm)
        #itp = LinearInterpolation(rowm, Ppj)
        #Pdk = itp(detm)
        
        while detm[detInd+c2]<rowm[end]
            sinoTmp[detInd] = sinoTmp[detInd]+(Pdk[detInd+c2]-Pdk[detInd])/deltaDetm
            detInd = detInd + detIinc
        end
        sinoTmp[detInd] = sinoTmp[detInd]+(Ppj[end]-Pdk[detInd])/deltaDetm                   
        
    end # Row loop 
    
    sinogram[proj,:] = sinoTmp .* L

end # Projection loop

return sinogram
end #endfunc


## Backprojection Branchless
function backprojection(sinogram,geo ; draw::Bool = false)

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
pixelX = pixelX .* pSize
pixelY = reverse(pixelY .- pSize/2, dims = 1)

# Tube
tubeX = 0     
tubeY = DSO    

# Iso-center
isoX = geo.isoX
isoY = geo.isoY

reconImg = zeros(nPix,nPix)
reconImgTmp = reconImg

# For each projection
for proj in 1:size(theta,1)
    
    angle = theta[proj]
       
    # Tubre rotation
    rtubeX = ( (tubeX - isoX )*cos(angle) - (tubeY - isoY )*sin(angle) ) + isoX
    rtubeY = ( (tubeX - isoX )*sin(angle) + (tubeY - isoY )*cos(angle) ) + isoX

    # Detector rotation
    rdetX = ( (detX .- isoX ).*cos(angle) - (detY .- isoY ).*sin(angle) ) .+ isoX
    rdetY = ( (detX .- isoX ).*sin(angle) + (detY .- isoY ).*cos(angle) ) .+ isoX
    
    if (draw == true)
        drawGeo(rtubeX,rtubeY,rdetX,rdetY)
    end    
    
    # Define angle case & which axis it it project boundaries
    # Case 1 
    if (((angle>=0)&&(angle<=pi/4))||(angle>=7*pi/4))       
        axisXCase = true   # Map on X axis()
        angleCase = 1
        c1=1
        c2=0
    else
        # Case 2
        if ((angle>pi/4)&&(angle<3*pi/4))
            axisXCase = false   # Map on Y axis()
            angleCase = 2
            c1=1
            c2=0
        else
            # Case 3
            if (((angle>=3*pi/4)&&(angle<=5*pi/4)))
                axisXCase = true   # Map on X axis()
                angleCase = 3
                c1=0
                c2=1
            else
            # Case 4
                axisXCase = false   # Map on Y axis()
                angleCase = 4
                c1=0
                c2=1
            end
        end
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
    L1 = zeros(1, nPix) 
    
    if (axisXCase) 
        # X-Ray pixel intersection calculation
        L = abs(pSize/cos(angle)) # This account for parallel-beam
        # Correction for fan-beam
        for n in 1:nPix
            L1[n]=sqrt((rtubeX-pixm[center_det])^2+(rtubeY-0)^2)/sqrt((rtubeX-pixm[n])^2+(rtubeY-0)^2)
        end
    else
        # X-Ray pixel intersection calculation
        L = abs(pSize/sin(angle)) # This account for parallel-beam
        # Correction for fan-beam
        for n in 1:nPix
            L1[n]=sqrt((rtubeX-0)^2+(rtubeY-pixm[center_det])^2)/sqrt((rtubeX-0)^2+(rtubeY-pixm[n])^2)
        end           
    end     
    L = L./L1
    
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
        
        reconTmp = zeros(1,nPix)
        
        rowm = pixm[row,:] # Get first mapped row from image.
        
        detInd = detIstart
        pixInd = pixIstart            
        
        # Find first detector overlap maped with pixel maped [Case 1]
        if (detm[detInd]-rowm[pixIstart]<-deltaDetm)
            while ((detm[detInd]-rowm[pixIstart]<-deltaDetm))            
                detInd = detInd + detIinc            
            end
        else
        # Find first pixel overlap maped with detector maped [Case 2]           
            if (detm[detIstart]-rowm[pixInd]>deltaPixm)
                while (detm[detIstart]-rowm[pixInd]>deltaPixm && pixInd<=nPix)            
                    pixInd = pixInd + pixIinc           
                end
            end
        end
        
        Ppj = integrate1D(sinogram[proj,:],deltaDetm) #Ppj is a row vector
        Pdk = Spline1D(sort!(detm), vec(Ppj); k=1)(rowm)
        #itp = LinearInterpolation(vec(detm), Ppj)
        #Pdk = itp(rowm)
        
        while (pixInd < nPix+1)
            reconTmp[pixInd] = reconTmp[pixInd]+(Pdk[pixInd+c1]-Pdk[pixInd+c2])/deltaDetm
            pixInd = pixInd + pixIinc 
        end
        # reconTmp[pixInd] = reconTmp[pixInd]+(Ppj[end]-Pdk[pixInd])/deltaPixm;   
        
        reconImgTmp[row,:] = reconTmp .* L
        
    end # Row loop 
    
    if ((angleCase == 1)||(angleCase == 3))
        reconImg = reconImg + reconImgTmp 
    else
        reconImg = reconImg + reverse(reconImgTmp, dims = 2)'
    end
   
   reconImg = reconImg ./ proj     
end # Projection loop

return reconImg
end #endfunc



## Integrate function
function integrate1D(p_v,pixelSize) # return Ppj

    n_pixel = size(p_v,1) #finds the number of rows in p_v
    
    P_x = 0
    Ppj = zeros(1,n_pixel+1)
    
    for pj in 1:n_pixel
        
       P_x = P_x + p_v[pj] * pixelSize
       
       Ppj[pj+1] = P_x
       
    end

    #Ppj[1] = -0.5 * P_x
    #Ppj = Ppj.-Ppj[1]
    return Ppj
end

## Map Y
# Function that map detector | pixel bounderies onto Y axis()
function mapp2y(x1,y1,x2,y2)
    y = y1-x1*(y1-y2)/(x1-x2)
    return y
end

## Map X
# Function that map detector | pixel bounderies onto X axis()
function mapp2x(x1,y1,x2,y2)
    x = -(y1)*(x1-x2)/(y1-y2)+x1
    return x
end

## Draw 2D
# Function to draw geometry
function drawGeo(tubeX,tubeY,detX,detY)
    plot(tubeX,tubeY,markershape = :asterisk)
    plot!([detX[1],detX[end]],[detY[1],detY[end]])
end


@show size(phantomImg)
jim(phantomImg)


# Projection
sinogramB = projectionBranchless(phantomImg,geo)
@show size(sinogramB)
p1 = jim(sinogramB)

#=
# Back-projection
imageB = backprojection(sinogramB,geo)
@show size(imageB)
@show imageB
p2 = jim(imageB)

jim(p1, p2)
=#




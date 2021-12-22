export fbp2_window

using FFTW

"""
    window = fbp2_window(n,window)


compute an apodizing window of length n and fft shift it
"""
function fbp2_window(n::Int,window::Symbol)

    # boxcar?
    if window===:none || window===:ramp
        window = ones(n)
        #=
    elseif window===:boxcar
        #cut = sscanf(window, 'boxcar,%g');
        window = my_boxcar(n, cut)

    elseif window===:hamming
        #cut = sscanf(window, 'hamming,%g');
        window = my_hamming(n, cut)

    elseif window===:hanning
        #cut = sscanf(window, 'hanning,%g');
        window = my_hann(n, cut)
        =#
    elseif window===:hann
        window = my_hann(n, 1.0)
    elseif window===:hann50
        window = my_hann(n, 0.5)
    elseif window===:hann75
        window = my_hann(n, 0.75)
    elseif window===:hann80
        window = my_hann(n, 0.80)

    else
        throw("unknown window $window")
    end

    return fftshift(window)
end

function window(n::Int,window::AbstractVector{<:Real})

    return fftshift(window)
end


function my_boxcar(n::Int, cutoff::Float64)
    w = round(cutoff * n)
    ii = (0:n-1) .- n/2
    return (abs.(ii) .< w/2)
    
end


function my_hann(n::Int, cutoff::Float64)
    w = round(cutoff * n)
    ii = (0:n-1) .- n/2
    return 0.5 * (1 .+ cos.(2*pi*ii/w)) .* (abs.(ii) .< w/2)
    
end


function my_hamming(n::Int, cutoff::Float64)
    w = round(cutoff * n)
    ii = (0:n-1) .- n/2
    return (0.54 .+ 0.46 * cos.(2*pi*ii/w)) .* (abs.(ii) .< w/2)
    
end



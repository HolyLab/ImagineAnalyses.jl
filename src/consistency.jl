using ImagineAnalyses
using Base.Test

#for measuring consistency of a periodic 1d timeseries
#returns a sample vector one cycle in duration where each sample is equal to the largest absolue difference measured between pairs of corresponding samples
#can choose to ignore the first ignore_ncycs cycles
function largest_cycle_diff(a::T, cycle_dur::Union{Int, HasTimeUnits}, delay_cyc1::Int=0; ignore_ncycs::Int=0) where {T<:AbstractVector}
    cyc_mat = get_cycles(a, cycle_dur, delay_cyc1)
    largest_cycle_diff(view(cyc_mat, :, ignore_ncycs+1:size(cyc_mat,2)))
end

largest_cycle_diff(a::T) where {T<:AbstractMatrix} = maximum(a, dims=2) .- minimum(a, dims=2)
largest_cycle_diff(sig::ImagineSignal, cycle_dur, delay_cyc1::Int; ignore_ncycs::Int=0) = largest_cycle_diff(get_samples(sig), cycle_dur, delay_cyc1; ignore_ncycs=ignore_ncycs)
largest_cycle_diff(sig::ImagineSignal, cycle_dur, delay_cyc1::HasTimeUnits; ignore_ncycs::Int=0) = largest_cycle_diff(get_samples(sig), cycle_dur, ImagineInterface.calc_num_samps(delay_cyc1, samprate(sig)); ignore_ncycs=ignore_ncycs)
largest_cycle_diff(sig::ImagineSignal, cycle_dur; ignore_ncycs::Int=0) = largest_cycle_diff(get_samples(sig), cycle_dur, 0; ignore_ncycs=ignore_ncycs)

#same as largest_cycle_diff but returns per-sample std
function largest_cycle_std(a::T, cycle_dur::Union{Int, HasTimeUnits}, delay_cyc1::Int; ignore_ncycs::Int=0) where {T<:AbstractVector}
    cyc_mat = get_cycles(a, cycle_dur, delay_cyc1)
    largest_cycle_std(view(cyc_mat, :, ignore_ncycs+1:size(cyc_mat,2)))
end

largest_cycle_std(a::T) where {T<:AbstractMatrix} = std(a, 2)
largest_cycle_std(sig::ImagineSignal, cycle_dur::HasTimeUnits, delay_cyc1::Int; ignore_ncycs::Int=0) = largest_cycle_std(get_samples(sig), cycle_dur, delay_cyc1; ignore_ncycs=ignore_ncycs)
largest_cycle_std(sig::ImagineSignal, cycle_dur::HasTimeUnits, delay_cyc1::HasTimeUnits; ignore_ncycs::Int=0) = largest_cycle_std(get_samples(sig), cycle_dur, ImagineInterface.calc_num_samps(delay_cyc1, samprate(sig)); ignore_ncycs=ignore_ncycs)
largest_cycle_std(sig::ImagineSignal, cycle_dur::HasTimeUnits; ignore_ncycs::Int=0) = largest_cycle_std(get_samples(sig), cycle_dur, 0; ignore_ncycs=ignore_ncycs)

#Reshape a periodic 1D signal into a matrix with sample index on the first axis and cycle index on the second
get_cycles(sig::ImagineSignal, cycle_dur, delay::HasTimeUnits) = get_cycles(get_samples(sig), cycle_dur, delay)
get_cycles(sig::ImagineSignal, cycle_dur, delay_nsamps::Int) = get_cycles(get_samples(sig), cycle_dur, delay_nsamps)
get_cycles(sig::ImagineSignal, cycle_dur) = get_cycles(get_samples(sig), cycle_dur, 0)

isapprox_int(x) = isapprox(round(Int, x), x)

get_cycles(a::AxisArray{T,1,D,Ax}, cycle_dur::HasTimeUnits; delay::HasTimeUnits = 0.0s) where {T,D,Ax} = get_cycles(a, cycle_dur, round(Int, delay/(step(axisvalues(a)[1]))))

function get_cycles(a::AxisArray{T,1,D,Ax}, cycle_dur::HasTimeUnits, delay_nsamps::Int) where {T,D,Ax}
    samp_dur = step(axisvalues(a)[1])
    nsamps = cycle_dur/samp_dur
    if isapprox_int(nsamps)
        nsamps = round(Int,nsamps) #special case
    end
    get_cycles(a.data, nsamps, delay_nsamps)
end

#general case (allocates new array)
#interpolates new cycle vector to honor the cycle time
#if input vector doesn't contain a whole number of cycles then omit a partial cycle at the end
function get_cycles(a::T, cycle_nsamps::Float64, delay_nsamps::Int = 0) where {T<:AbstractVector}
    la = length(a) - delay_nsamps
    if la < cycle_nsamps
        error("Insufficient sample count for the specified cycle time")
    end
    ncycs = la/cycle_nsamps
    if isapprox_int(ncycs)
        ncycs = round(Int, ncycs)
    else
        ncycs = floor(Int, ncycs)
    end
    nsamps_interp = ceil(Int, cycle_nsamps)
    output = zeros(eltype(a), nsamps_interp, ncycs)
    a_itp = extrapolate(interpolate(a, BSpline(Linear())), Flat())
    step_size = cycle_nsamps/nsamps_interp
    for i = 1:ncycs
        #output[:,i] = a_itp[linspace(1+(i-1)*cycle_nsamps, i*cycle_nsamps, nsamps_interp)] #not sure why this fails for certain float ranges
        cyc_inds = range(delay_nsamps+1+1+(i-1)*cycle_nsamps, stop=i*cycle_nsamps, length=nsamps_interp)
        for (ii,ival) in enumerate(cyc_inds)
            output[ii,i] = a_itp(ival)
        end
    end
    return output
end

#special case (in-place reshape)
function get_cycles(a::T, cycle_nsamps::Int, delay_nsamps::Int = 0) where {T<:AbstractVector}
    la = length(a) - delay_nsamps
    d, r = divrem(la, cycle_nsamps)
    if la < cycle_nsamps
        error("Insufficient sample count for the specified cycle time")
    else
        a = view(a, (delay_nsamps+1):(delay_nsamps+d*cycle_nsamps))
    end
    return reshape(a, cycle_nsamps, d)
end

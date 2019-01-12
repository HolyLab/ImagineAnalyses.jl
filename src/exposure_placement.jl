#Iterate through samps and circ list, recording the indices at which the samp vector crosses the thresholds in circ_list
#circ_list is expected to be sorted in either forward or reverse order, and the first index of samps is expected to be outside the range spanned by first(circ_list) and last(circ_list)
#After fully iterating through circ_list the search will continue on the reversed circ_list.
#Since the endpoints of circ_list get repeated they can be counted too often due to noise.  For this reason the pad_samps argument determines how long to wait before counting another threshold crossing each time the list is reversed
#Returns a vector-of-vectors where each vector contains the crossing indices of a particular circ_list value
function find_circular(samps::AbstractVector, circ_list::AbstractVector, pad_nsamps::Int)
    is_incr = true
    reversed = false
    rlist = reverse(circ_list)
    if !issorted(circ_list)
        if !issorted(rlist)
            error("circ_list must be in increasing or decreasing order")
        end
        circ_list = rlist
        is_incr = false
        reversed = true
    end
    samp1 = first(samps)
    if (is_incr && (samp1 > circ_list[1])) || (!is_incr && (samp1 < circ_list[end]))
        error("The sample vector begins in the middle of a cycle")
    end
    outputs = Vector{Vector{Int}}(undef, length(circ_list))
    for i = 1:length(circ_list)
        outputs[i] = Int[]
    end
    wait_stop = 0
    circ_i = reversed ? length(circ_list) : 1
    for i = 1:length(samps)
        if i == wait_stop
            if (is_incr && (samps[i] >= circ_list[circ_i])) || (!is_incr && (samps[i] <= circ_list[1]))
                error("The signal had already passed threshold after waiting for pad_nsamps samples")
            end
        elseif i > wait_stop
            if is_incr && (samps[i] >= circ_list[circ_i])
                push!(outputs[circ_i], i)
                if circ_i == length(circ_list)
                    is_incr = false
                    wait_stop = i + pad_nsamps
                else
                    circ_i += 1
                end
            elseif !is_incr && (samps[i] <= circ_list[circ_i])
                push!(outputs[circ_i], i)
                if circ_i == 1
                    is_incr = true
                    wait_stop = i + pad_nsamps
                else
                    circ_i -= 1
                end
            end
        end
    end
    if reversed
        outputs = reverse(outputs)
    end
    return outputs
end

#Returns a vector that is the mean of each cycle in "samps" beginning at start_idx
#Setting start_idx greater than one allows skipping of early samples while the piezo warms up
function mean_cycle(samps::AbstractVector, nsamps_cycle::Int; start_idx = 1)
    lsamps = length(samps)
    @assert start_idx >=1 && start_idx <= lsamps
    @assert (lsamps-start_idx+1) >= nsamps_cycle
    rng_analyzed = start_idx:lsamps
    samps_used = view(samps, rng_analyzed)
    if mod(length(samps_used), nsamps_cycle) != 0
        @warn "Found a fractional number of cycles.  Truncating to analyze a whole number of cycles."
        rng_analyzed = start_idx:start_idx+(nsamps_cycle*div(length(rng_analyzed),nsamps_cycle)-1)
        samps_used = view(samps, rng_analyzed)
    end
    cyc_mat = get_cycles(samps_used, nsamps_cycle)
    print("Averaging $(size(cyc_mat,2)) cycles...\n")
    return reshape(mean(ustrip.(cyc_mat), dims=2), nsamps_cycle)  * unit(samps[1]) #work around limitation of mean function with Unitful Quantities
end

#find integer sample shift to align mod and mon.  (Even with calibration there is often a phase delay between the MOD and MON signal)
#We do this so that we know when the MON signal reaches its minimum
#indmax(xcorr(a,b)) is displaced from the center of the xcorr vector by the amount that b would need to be shifted to align with a
function mon_delay(mod_cyc, mon_cyc)
    xc = xcorr(repeat(ustrip.(mon_cyc), 3), ustrip.(mod_cyc)) #replicate it so that edges don't corrupt results
    ctr_i = div(length(xc),2)+1
    half_window_sz = div(length(mon_cyc),4)
    inner_half = xc[(ctr_i-half_window_sz):(ctr_i+half_window_sz)]
    #amount that mod_cyc needs to be shifted to align with mon_cyc (should be positive)
    #we examine only the inner half of shifts because we don't want to shift more than a half cycle
    return argmax(inner_half) - (div(length(inner_half),2) + 1)
end

#Returns a set of indices marking when the mon signal crosses the values in the crossings vector.
#Useful for deciding when to time camera exposures / laser pulses
#note: if calibration is needed it must be done before calling this function. i.e. sig should be calibrated already
function pulse_timings(mod_cyc::AbstractVector, mon_cyc::AbstractVector, targets::AbstractVector, pad_nsamps::Int)
    @assert length(mod_cyc) == length(mon_cyc)
    nsamps_cycle = length(mod_cyc)
    mon_idxs = find_circular(mon_cyc, targets, pad_nsamps)
    @assert all(map(length, mon_idxs).==2) #should have exactly two crossings for each target with a cyclical waveform
    return mon_idxs
end

function flash_cam_cycs(mod_cyc::AbstractVector, mon_cyc::AbstractVector, slice_zs, pad_nsamps::Int, flash_time::HasTimeUnits, exp_time::HasTimeUnits, sample_rate::HasInverseTimeUnits; is_bidi = true)
    offset = mon_delay(ustrip.(mod_cyc), ustrip.(mon_cyc))
    mon_cyc = circshift(mon_cyc, -offset)
    @assert flash_time <= exp_time
    flash_ctr_is = pulse_timings(mod_cyc, mon_cyc, slice_zs, pad_nsamps)
    flash_ctr_fwd = [x[1] for x in flash_ctr_is]
    flash_ctr_back = reverse([x[2] for x in flash_ctr_is]) #reversed so that corresponding z locations line up
    if is_bidi
        flash_ctr_is = vcat(flash_ctr_fwd, flash_ctr_back)
    else
        flash_ctr_is = flash_ctr_fwd
    end
    flash_cam_cycs(flash_ctr_is, flash_time, exp_time, sample_rate, length(mon_cyc), offset)
end

function flash_cam_cycs(flash_ctr_is, flash_time::HasTimeUnits, exp_time::HasTimeUnits, sample_rate::HasInverseTimeUnits, nsamps_cyc::Int, offset::Int)
    #@assert all(map((x,y) -> isapprox(x, y, rtol = (1/typemax(Int16) * 800μm)), (mean_cyc[flash_ctr_fwd], mean_cyc[flash_ctr_back])))
    min_flash_sep = minimum(diff(flash_ctr_is))
    max_exp_time = ((min_flash_sep-1) / sample_rate) * 0.97
    #@show max_exp_time = (1/stack_rate - (1/samprate(pos) *length(slice_zs))) / length(slice_zs) #leaves one low sample between exposure pulses
    if exp_time >= max_exp_time
        @warn "Could not use the requested exposure time while keeping sufficient separation between exposures.  Using an exposure time of $(max_exp_time) instead.  This means that the maximum framerate required from the camera will be $(100.0 * inv(max_exp_time)/inv(exp_time))% of the requested rate. (choose the ROI size to meet this requirement)"
        exp_time = max_exp_time
    end
    #now place flashes, centered on the indices above
    nsamps_flash = ImagineInterface.calc_num_samps(flash_time, sample_rate)
    if iseven(nsamps_flash)
	    nsamps_flash -= 1 #force odd so that fwd and reverse flashes line up in bidi recordings
    end
    nsamps_from_ctr = div(nsamps_flash,2)
    flash_ivs = map(x-> ClosedInterval(x-nsamps_from_ctr, x+nsamps_from_ctr), flash_ctr_is)
    #now place exposures, aligning exposure ends with flash ends
    nsamps_exp = ImagineInterface.calc_num_samps(exp_time, sample_rate)
    exp_ivs= map(x-> ClosedInterval(x+nsamps_from_ctr-nsamps_exp+1, x+nsamps_from_ctr), flash_ctr_is)
    flash_cyc = ImagineInterface.gen_pulses(nsamps_cyc, flash_ivs)
    cam_cyc = ImagineInterface.gen_pulses(nsamps_cyc, exp_ivs)

    return flash_cyc, cam_cyc, offset
end

#empirical max acceleration
max_acc(v, sr) = maximum(abs.(upreferred.(diff(diff(v)./(1/sr))./(1/sr))))
max_vel(v, sr) = maximum(abs.(upreferred.(diff(v)./(1/sr))))

#Returns a triangle wave equivalent to the input (same max, min, freq), useful for comparisons
function eq_triangle(v, sr::HasInverseTimeUnits)
    mx = maximum(v)
    mn = minimum(v)
    freq = upreferred(sr/length(v))
    fwd, bck = ImagineInterface.gen_bidi_pos(mn, mx, upreferred(1/(2.0*freq)), sr)
    return vcat(fwd,bck)
end

#calculate the fraction of the camera's max framerate that can be utilized with the given sample vector
#by comparing the maximum velocity reached with the ideal constant velocity of a triangle wave
function framerate_efficiency(v, sr::HasInverseTimeUnits)
    tri = eq_triangle(v, sr)
    vin = max_vel(v,sr)
    vtri = max_vel(tri,sr)
    return upreferred(vtri/vin)
end

#calculate the fraction of the camera's max framerate that can be utilized with the given sample vector and slice locations
#This version is more practical than above (above is more of a worst-case efficiency)
function framerate_efficiency(v, sr::HasInverseTimeUnits, slice_zs; pad_nsamps=ImagineInterface.calc_num_samps(0.0005s, sr), toffsets_fwd=fill(0.0s, length(slice_zs)), toffsets_bck=fill(0.0s, length(slice_zs)))
    mexp = max_exp(v,sr,slice_zs; pad_nsamps=pad_nsamps, toffsets_fwd=toffsets_fwd, toffsets_bck=toffsets_bck)
    tri = eq_triangle(v, sr)
    mexp_tri = max_exp(tri, sr, slice_zs; pad_nsamps=2)
    @assert eltype(v) == eltype(slice_zs)
    return upreferred(mexp/mexp_tri)
end

#calculate longest exposure time possible given sample vector v
#Also considers temporal offsets that have been found empirically (toffsets_fwd and toffsets_bck keyword args)
function max_exp(v, sr::HasInverseTimeUnits, slice_zs; pad_nsamps = ImagineInterface.calc_num_samps(0.0005s, sr), toffsets_fwd=fill(0.0s, length(slice_zs)), toffsets_bck = fill(0.0s,length(slice_zs)))
    first_idx = argmin(v)
    v = circshift(v, -first_idx) #assumes slice_zs are sorted and increasing
    timings = ImagineAnalyses.find_circular(v, slice_zs, pad_nsamps)
    fwdt = [x[1] for x in timings]
    bckt = reverse([x[2] for x in timings])
    toffsets_fwd_nsamps = map(x->ImagineInterface.calc_num_samps(x, sr), toffsets_fwd)
    toffsets_bck_nsamps = map(x->ImagineInterface.calc_num_samps(x, sr), toffsets_bck)
    fwdt = fwdt.-toffsets_fwd_nsamps
    bckt = bckt.-reverse(toffsets_bck_nsamps)
    diffi = minimum(diff(vcat(fwdt, bckt))) - 1 #subtract one because we need a sample between exposures (to send TTL down)
    diffi = min(diffi, first(fwdt)+(length(v)-last(bckt)) - 1) #handle first/last timing
    return upreferred(diffi / sr)
end

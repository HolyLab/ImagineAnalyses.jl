#__precompile__(true)

module ImagineAnalyses

using AxisArrays, Interpolations, Unitful, IntervalSets
using Reexport
@reexport using CachedCalls, ImagineInterface, UnitAliases

import UnitAliases.HasTimeUnits
import CachedCalls.CachedCall
import Unitful: μm, s, Hz
using DSP
using Statistics

export dumb_bench,
        get_cycles,
        largest_cycle_diff,
        largest_cycle_std,
        mean_cycle,

        find_circular,
        get_circular,
        mon_delay,
        pulse_timings,
        flash_cam_cycs,
        max_exp,
        framerate_efficiency


#This will be a bit of a hodgepodge of analyses one can do on ImagineSignals and acquired images.  May be worth organizing differently if it gets too big.

#this may do something like measure the lag between sig1 and sig2
dumb_bench(sig1::ImagineSignal, sig2::ImagineSignal) = 0.0

#Define new constructors for CachedCalls on a per-function basis when the return type is known (warning: risky when function isn't type stable)
CachedCall(f::typeof(dumb_bench), args...) = CachedCall("dumb benchmark", f, (args...,), Float64)

include("consistency.jl")
include("exposure_placement.jl")

#To implement:

#lag(sig1::ImagineSignal, sig2::ImagineSignal) 
    #lag(sig1::AxisArray, sig2::AxisArray) 
#offset(sig1::ImagineSignal, sig2::ImagineSignal) 
#gain(sig_out::ImagineSignal, sig_in::ImagineSignal)
#corr(sig1::ImagineSignal, sig2::ImagineSignal)
#pulse_count(sig::ImagineSignal)
#range(sig::ImagineSignal)
#signal_consistency(sig::ImagineSignal, cycle_dur::HasTimeUnits) #consistency of a periodic 1d timeseries. may change the interface

#The below use images
#frame_count(camsig::ImagineSignal, img::AbstractImage{T,N}) where {T,N}
#frame_order(camsig::ImagineSignal, img::AbstractImage{T,N}) where {T,N}  #returns an integer sequence of frame numbers
#bidi_image_consistency(img::AbstractImage{T,4}) where {T} #returns each pair of ssds in a bidi sequence
#image_consistency(img::AbstractImage{T,N}, imgtemplate::AbstractImage{T,N-1}) where {T, N} #returns ssd between imgtemplate and each image in the img sequence

end # module

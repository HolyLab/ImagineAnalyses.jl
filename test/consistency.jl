using ImagineAnalyses
using Test

using ImagineInterface, Unitful, Interpolations, AxisArrays
import Unitful:s, μm

###get_cycles

##whole number of cycles
#raw array
cyc = [1.0:1.0:10.0...]
ncycs = 5
cyc_mat = repeat(cyc, 1, ncycs)
cyc_vec = reshape(cyc_mat, (ncycs * length(cyc),))
cyc_mat2 = ImagineAnalyses.get_cycles(cyc_vec, length(cyc))
@test cyc_mat == cyc_mat2 #test special case with whole number of cycles
#AxisArray
#(with sample-aligned cycles)
cyc_vec_aa = AxisArray(cyc_vec, Axis{:time}(range(1.0s, stop=ncycs*10.0s, length=ncycs*length(cyc))))
@test all(ImagineAnalyses.get_cycles(cyc_vec_aa, 10.0s) .== cyc_mat)

#(without sample-aligned cycles)
#cycle still takes 10 seconds, but we sample at 1.22 samples per second which works out evenly to 5 cycles in 61 samples
cyc_vec_itp = extrapolate(interpolate(cyc_vec, BSpline(Linear())), Flat())
#cyc_vec3 = cyc_vec_itp(range(1, stop=length(cyc_vec), length=61))
cyc_vec3 = similar(cyc_vec, (61,))
#for (i, t) in enumerate(range((1/1.22), stop=length(cyc_vec), length=61))
#cyc_vec3 = cyc_vec_itp(range(1/1.22, stop=length(cyc_vec), length=61)) #not sure why this doesn't work
for (i, t) in enumerate(range(1/1.22, stop=length(cyc_vec), length=61))
    cyc_vec3[i] = cyc_vec_itp(t)
end
cyc_vec3_aa = AxisArray(cyc_vec3, Axis{:time}(range((1.0s/1.22), stop=50.0s, length=length(cyc_vec3))))
cycs_out = ImagineAnalyses.get_cycles(cyc_vec3_aa, 10.0s)
@test size(cycs_out) == (13,5) #12.2 samples per cycle but it rounds up
same_subcyc = cycs_out[3:11,:] #due to interpolation at edges only the middle portion of the cycle vector is consistent
for i = 1:size(same_subcyc,1)
    val = same_subcyc[i,1]
    @test all(isapprox.(val,same_subcyc[i,:]))
end

##partial cycle
#raw array
cyc = [1.0:1.0:10.0...]
ncycs = 5
cyc_mat = repeat(cyc, 1, ncycs)
cyc_vec = reshape(cyc_mat, (ncycs * length(cyc),))
cyc_vec2 = push!(copy(cyc_vec), 1.0)
@test all(ImagineAnalyses.get_cycles(cyc_vec2, length(cyc)) .== cyc_mat) #should ignore the last partial cycle
#AxisArray
#(with sample-aligned cycles)
cyc_vec_aa = AxisArray(push!(cyc_vec,1.0), Axis{:time}(range(1.0s, stop=51.0s, length=51)))
@test all(ImagineAnalyses.get_cycles(cyc_vec_aa, 10.0s) .== cyc_mat)

#(without sample-aligned cycles)

#cycle still takes 10 seconds, but we sample at 1.02 samples per second which does not work out evenly
cyc_vec_itp = extrapolate(interpolate(cyc_vec, BSpline(Linear())), Flat())
#cyc_vec3 = cyc_vec_itp(range(1, stop=length(cyc_vec), length=61))
cyc_vec3 = similar(cyc_vec, (51,))
#for (i, t) in enumerate(range((1/1.22), length(cyc_vec), 61))
for (i, t) in enumerate(range(1/1.02, stop=51*(1/1.02), length=51))
    cyc_vec3[i] = cyc_vec_itp(t)
end
cyc_vec3_aa = AxisArray(cyc_vec3, Axis{:time}(range((1.0s/1.02), stop=51.0*(1s/1.02), length=length(cyc_vec3))))
cycs_out = ImagineAnalyses.get_cycles(cyc_vec3_aa, 10.0s)
@test size(cycs_out) == (11,5)
same_subcyc = cycs_out[2:10,:] #due to interpolation at edges only the middle portion of the cycle vector is consistent
for i = 1:size(same_subcyc,1)
    val = same_subcyc[i,1]
    @test all(isapprox.(val,same_subcyc[i,:]))
end

#ImagineSignal
ocpi2 = rigtemplate("ocpi-2"; sample_rate=10200*inv(Unitful.s)) #1.02 * 10k
pos = first(getpositioners(ocpi2))
append!(pos, "test", cyc_vec3*Unitful.μm)
samps_vec = get_samples(pos) #returns an AxisArray
sig_cycs = ImagineAnalyses.get_cycles(pos, (10.0s/10000))
@test all(sig_cycs.==ImagineAnalyses.get_cycles(samps_vec, 10.0s/10000))

#largest_cycle_diff
cyc = [1.0:1.0:10.0...]
ncycs = 5
cyc_mat1 = repeat(cyc, 1, ncycs)
cyc_mat2 = repeat(cyc, 1, ncycs)
cyc_mat2[:,3] .+= 1.0

cyc_vec1 = reshape(cyc_mat1, (ncycs * length(cyc),))
cyc_vec2 = reshape(cyc_mat2, (ncycs * length(cyc),))
@test all(largest_cycle_diff(cyc_vec1, 10).==0.0)
@test all(largest_cycle_diff(cyc_vec2, 10).==1.0)
@test all(largest_cycle_diff(cyc_vec2, 10; ignore_ncycs=3).==0.0) #ignore first 3

#AxisArray
cyc_vec2_aa = AxisArray(cyc_vec2, Axis{:time}(range(1.0s, stop=ncycs*10.0s, length=ncycs*length(cyc))))
@test all(largest_cycle_diff(cyc_vec2_aa, 10.0s).==1.0)
@test all(largest_cycle_diff(cyc_vec2_aa, 10.0s; ignore_ncycs=3).==0.0) #ignore first 3

#ImagineSignal
ocpi2 = rigtemplate("ocpi-2"; sample_rate=1*inv(Unitful.s))
pos = first(getpositioners(ocpi2))
append!(pos, "test", cyc_vec2*Unitful.μm)
myapprox = (x,y)->isapprox(ustrip(x),ustrip(y); rtol=sqrt(1/typemax(Int16))) #account for Int16-limited encoding of piezo signal
@test all(map(x->myapprox(x,1.0μm),  largest_cycle_diff(pos, 10.0s)))
@test all(map(x->myapprox(x,0.0μm),  largest_cycle_diff(pos, 10.0s; ignore_ncycs=3)))

samps_vec = get_samples(pos)
ocpi2 = rigtemplate("ocpi-2"; sample_rate=1*inv(Unitful.s))
pos = first(getpositioners(ocpi2))
samps_vec = vcat(fill(0.0*Unitful.μm, 111), samps_vec)
append!(pos, "test", samps_vec)
@test all(map(x->myapprox(x,1.0μm),  largest_cycle_diff(pos, 10.0s, 111))) #111 sample delay period
@test all(map(x->myapprox(x,1.0μm),  largest_cycle_diff(pos, 10.0s, 111/samprate(pos)))) #time units


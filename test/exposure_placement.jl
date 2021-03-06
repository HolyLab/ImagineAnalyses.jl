using ImagineAnalyses, IntervalSets, Unitful, Test
import Unitful:μm, s, Hz

#find_circular
samps = repeat([1:9...; 10:-1:2...], 3)
circ_list = [2; 4]
pad_nsamps = 0
i_s = find_circular(samps, circ_list, pad_nsamps)
@test i_s[1] == findall(x->x==circ_list[1], samps)
@test i_s[2] == findall(x->x==circ_list[2], samps)
@test all(find_circular(samps, circ_list, 1) .== i_s)
@test_throws Exception find_circular(samps, circ_list, 2)
circ_list = [4; 2]
@test_throws Exception find_circular(samps, circ_list, pad_nsamps) #because it doesn't begin decreasing
samps = repeat([10:-1:2...; 1:9...], 3) #begins decreasing
i_s = find_circular(samps, circ_list, pad_nsamps)
@test i_s[1] == findall(x->x==circ_list[1], samps)
@test i_s[2] == findall(x->x==circ_list[2], samps)

samps = repeat([1:9...; 10:-1:2...], 3)
circ_list = [2.5; 4.0] #fractional
i_s = find_circular(samps, circ_list, pad_nsamps)
@test i_s[1][1:2:end] == findall(x->x==3, samps)[1:2:end] #on increasing sweep
@test i_s[1][2:2:end] == findall(x->x==2, samps)[2:2:end] #on decreasing sweep
@test i_s[2] == findall(x->x==circ_list[2], samps) #unaffected

#again with units
samps *= Unitful.μm
circ_list *= Unitful.μm
i_s = find_circular(samps, circ_list, pad_nsamps)
@test i_s[1][1:2:end] == findall(x->x==3μm, samps)[1:2:end] #on increasing sweep
@test i_s[1][2:2:end] == findall(x->x==2μm, samps)[2:2:end] #on decreasing sweep
@test i_s[2] == findall(x->x==circ_list[2], samps) #unaffected

samps = fill(1.0μm, 30)
@test all(mean_cycle(samps, 10; start_idx = 1).==samps[1:10])
samps[1] = 2.0μm
@test all(mean_cycle(samps, 10; start_idx = 2).==samps[2:11])
samps[1:2] .= 2.0μm
@test all(mean_cycle(samps, 8; start_idx = 3).==samps[3:10])

pulse_mod = Float64.(ImagineInterface.gen_pulses(20, [10..12]))
pulse_mon = Float64.(ImagineInterface.gen_pulses(20, [12..14]))
@test mon_delay(pulse_mod, pulse_mon) == 2
@test mon_delay(pulse_mon, pulse_mod) == -2
@test mon_delay(convert(Vector{Float64}, pulse_mod).*Unitful.μm, convert(Vector{Float64}, pulse_mon).*Unitful.μm) == 2

pos_fwd, pos_back = ImagineInterface.gen_bidi_pos(0.0μm, 20.0μm, 1.0s, 20.0Hz)
mod_cyc = ustrip.(vcat(pos_fwd, pos_back))
mod_cyc = round.(Int, mod_cyc) * unit(pos_fwd[1]) #inexact nature of floats can cause make tests nonintuitive
mon_cyc = mod_cyc
targets = [5μm; 10μm; 15μm]
pad_nsamps = 0
pt = pulse_timings(mod_cyc, mon_cyc, targets, pad_nsamps)
for i=1:length(targets)
    @test all(mon_cyc[pt[i]].==targets[i])
end

flash_cyc, cam_cyc, offset = flash_cam_cycs(mod_cyc, mon_cyc, targets, pad_nsamps, 0.05s, 0.1s, 20.0Hz; is_bidi = true)
@test length(findall(flash_cyc)) == sum([length(x) for x in pt])
@test length(findall(cam_cyc)) == sum([length(x) for x in pt]) * 2
for i=1:length(targets)
    for ii in pt[i]
        @test flash_cyc[ii] == true
        @test all(cam_cyc[ii-1:ii].== true) #Exposure lasts two samples, begins the sample previous to the flash
    end
end

#now with a shift
mon_cyc = circshift(mod_cyc, 1)
pt = pulse_timings(mod_cyc, mon_cyc, targets, pad_nsamps)
for i=1:length(targets)
    @test all(mon_cyc[pt[i]].==targets[i])
end
flash_cyc, cam_cyc, offset = flash_cam_cycs(mod_cyc, mon_cyc, targets, pad_nsamps, 0.05s, 0.1s, 20.0Hz; is_bidi = true)
@test offset == 1
flash_cyc = circshift(flash_cyc, offset)
cam_cyc = circshift(cam_cyc, offset)
@test length(findall(flash_cyc)) == sum([length(x) for x in pt])
@test length(findall(cam_cyc)) == sum([length(x) for x in pt]) * 2
for i=1:length(targets)
    for ii in pt[i]
        @test flash_cyc[ii] == true
        @test all(cam_cyc[ii-1:ii].== true) #Exposure lasts two samples, begins the sample previous to the flash
    end
end

#function cam_flash_cycs(flash_ctr_is, flash_time::HasTimeUnits, exp_time::HasTimeUnits, sample_rate::HasInverseTimeUnits)

v = [1.0:1.0:10.0...;9.0:-1.0:2.0...].*Unitful.μm
zs = [2.0; 4.0; 6.0] .* Unitful.μm
sr = 1.0Hz

@test max_exp(v, sr, zs) == 1.0s
eff = framerate_efficiency(v, sr, zs)
@test eff == 1.0
@test framerate_efficiency(v, sr) == 1.0 #version that doesn't consider slice locations (worst case)

zs = [2.0; 5.0] .* Unitful.μm #last slice is limiting
@test max_exp(v, sr, zs) == 1.0s
eff = framerate_efficiency(v, sr, zs)
@test eff == 1.0

v = [1.0;3.0;5.0;6.0:10.0...;9.0:-1.0:2.0...].*Unitful.μm #nonlinear in the beginning of the forward sweep
zs = [3.0; 6.0] .* Unitful.μm #max exp should be 2.0 in the linear case, 1.0 due to nonlinearity: efficiency of 0.5
@test max_exp(v, sr, zs) == 1.0s
eff = framerate_efficiency(v, sr, zs)
@test eff == 0.5


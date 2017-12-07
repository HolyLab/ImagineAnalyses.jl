using ImagineAnalyses, Base.Test
import Unitful:μm, s

#find_circular
samps = repmat([1:9...; 10:-1:2...], 3)
circ_list = [2; 4]
pad_nsamps = 0
i_s = find_circular(samps, circ_list, pad_nsamps)
@test i_s[1] == find(x->x==circ_list[1], samps)
@test i_s[2] == find(x->x==circ_list[2], samps)
@test all(find_circular(samps, circ_list, 1) .== i_s)
@test_throws Exception find_circular(samps, circ_list, 2)
circ_list = [4; 2]
@test_throws Exception find_circular(samps, circ_list, pad_nsamps) #because it doesn't begin decreasing
samps = repmat([10:-1:2...; 1:9...], 3) #begins decreasing
i_s = find_circular(samps, circ_list, pad_nsamps)
@test i_s[1] == find(x->x==circ_list[1], samps)
@test i_s[2] == find(x->x==circ_list[2], samps)

samps = repmat([1:9...; 10:-1:2...], 3)
circ_list = [2.5; 4.0] #fractional
i_s = find_circular(samps, circ_list, pad_nsamps)
@test i_s[1][1:2:end] == find(x->x==3, samps)[1:2:end] #on increasing sweep
@test i_s[1][2:2:end] == find(x->x==2, samps)[2:2:end] #on decreasing sweep
@test i_s[2] == find(x->x==circ_list[2], samps) #unaffected

#again with units
samps *= Unitful.μm
circ_list *= Unitful.μm
i_s = find_circular(samps, circ_list, pad_nsamps)
@test i_s[1][1:2:end] == find(x->x==3μm, samps)[1:2:end] #on increasing sweep
@test i_s[1][2:2:end] == find(x->x==2μm, samps)[2:2:end] #on decreasing sweep
@test i_s[2] == find(x->x==circ_list[2], samps) #unaffected

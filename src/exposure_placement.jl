#Iterate through samps and circ list, marking the indices at which samp the samp vector crosses th thresholds in circ_list
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
    outputs = Vector{Vector{Int}}(length(circ_list))
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
    return reversed ? reverse(outputs) : outputs
end

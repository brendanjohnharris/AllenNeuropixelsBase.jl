struct HybridSession <: AbstractNWBSession
    file
    pyObject
end










"""
Get the lfp data for a probe, providing *indices* for channels and times. See function below for indexing by channel ids and time values/intervals
"""
function _getlfp(session::AbstractSession, probeid::Int; channelidxs=1:length(getlfpchannels(session, probeid)), timeidxs=1:length(getlfptimes(session, probeid)))
    @assert(any(getprobeids(session) .== probeid), "Probe $probeid does not belong to session $(getid(session))")
    @assert(subset(getprobes(session), :id=>ByRow(==(probeid)))[!, :has_lfp_data][1], @error "Probe $probeid does not have LFP data")
    path = getlfppath(session, probeid)
    if !isfile(path)
        downloadlfp(session, probeid)
    end
    if !isfile(path)
        @error "LFP data did not download correctly"
    end

    f = h5open(path)

    timedata = getlfptimes(session, probeid, timeidxs)
    dopermute = true
    channelids = getlfpchannels(session, probeid)
    channelids = channelids[channelidxs]
    if (channelidxs isa Union{Int64, AbstractRange{Int64}}) & (timeidxs isa Union{Int64, AbstractRange{Int64}}) # Can use HDF5 slicing
        lfp = f["acquisition"][splitext(basename(path))[1]][splitext(basename(path))[1]*"_data"]["data"][channelidxs, timeidxs]
    elseif timeidxs isa Union{Int64, AbstractRange{Int64}}
        lfp = [f["acquisition"][splitext(basename(path))[1]][splitext(basename(path))[1]*"_data"]["data"][i, timeidxs] for i ∈ channelidxs]
        lfp = hcat(lfp...)
        dopermute = false
    else
        lfp = read(f["acquisition"][splitext(basename(path))[1]][splitext(basename(path))[1]*"_data"]["data"])
        lfp = lfp[channelidxs, timeidxs]
    end
    if lfp isa Vector
       lfp = reshape(lfp, 1, length(lfp))
    end
    if channelids isa Number
        channelids = [channelids]
    end
    if dopermute
        lfp = permutedims(lfp, reverse(1:ndims(lfp)))
    end
    X = DimArray(lfp, (Ti(timedata),  Dim{:channel}(channelids)))
    close(f)
    return X
end

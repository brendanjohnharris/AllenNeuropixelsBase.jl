export HybridSession
using NWBStream
using PythonCall
using NWBStream
using DataFrames
using Downloads
using JSON
import NWBStream.url2df

mutable struct HybridSession <: AbstractNWBSession
    url::Any
    file::Any
    io::Any
    pyObject::Any
    function HybridSession(url::String, file, io, pyObject = ())
        S = new(url, file, io, pyObject)
        f(S::HybridSession) = @async s3close(S.io)
        finalizer(f, S)
    end
end
HybridSession(url::String) = (S = HybridSession(url, s3open(url)...); initialize!(S); S)
function HybridSession(session_id::Int, args...; kwargs...)
    HybridSession(VisualBehavior.getsessionfile(session_id), args...; kwargs...)
end

# Calls to getlfp should use the streaming approach if the length is less than some amount, otherwise download the file using the allensdk
function _loadlfp(session::HybridSession, probeid::Int;
                  channelidxs = 1:length(getlfpchannels(session, probeid)),
                  timeidxs = 1:length(getlfptimes(session, probeid)))
    @assert(any(getprobeids(session) .== probeid),
            "Probe $probeid does not belong to session $(getid(session))")
    @assert(subset(getprobes(session), :id => ByRow(==(probeid)))[!, :has_lfp_data][1],
            @error "Probe $probeid does not have LFP data")
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
    if (channelidxs isa Union{Int64, AbstractRange{Int64}}) &
       (timeidxs isa Union{Int64, AbstractRange{Int64}}) # Can use HDF5 slicing
        lfp = f["acquisition"][splitext(basename(path))[1]][splitext(basename(path))[1] * "_data"]["data"][channelidxs,
                                                                                                           timeidxs]
    elseif timeidxs isa Union{Int64, AbstractRange{Int64}}
        lfp = [f["acquisition"][splitext(basename(path))[1]][splitext(basename(path))[1] * "_data"]["data"][i,
                                                                                                            timeidxs]
               for i in channelidxs]
        lfp = hcat(lfp...)
        dopermute = false
    else
        lfp = read(f["acquisition"][splitext(basename(path))[1]][splitext(basename(path))[1] * "_data"]["data"])
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
    X = ToolsArray(lfp, (𝑡(timedata), Chan(channelids));
                   metadata = Dict(:sessionid => getid(session), :probeid => probeid))
    close(f)
    return X
end

function _streamlfp(session::AbstractNWBSession, probeid::Int;
                    channelidxs = 1:length(getlfpchannels(session, probeid)),
                    timeidxs = 1:length(getlfptimes(session, probeid)))
    @assert(any(getprobeids(session) .== probeid),
            "Probe $probeid does not belong to session $(getid(session))")
    _timeidxs = timeidxs .- 1 # Python sucks
    _channelidxs = channelidxs .- 1 # Python sucks
    f, io = getprobefile(session, probeid) |> NWBStream.s3open
    _lfp = Dict(f.acquisition)["probe_1108501239_lfp_data"]

    # timedata = getlfptimes(session, probeid)
    timedata = getlfptimes(session, probeid)
    timedata = timedata[timeidxs]
    channelids = getlfpchannels(session, probeid)
    channelids = channelids[channelidxs]
    lfp = zeros(Float32, length(timedata), length(channelids))
    for (i, _channelidx) in enumerate(_channelidxs)
        d = diff(_timeidxs)
        # if all(d == d[1]) # We can do some efficient slicing
        #     slc = @py slice(first(_timeidxs), last(_timeidxs), d[1])
        #     @time @py _lfp.data[slc]
        # else
        lfp[:, i] .= pyconvert(Vector{Float32}, _lfp.data[_timeidxs, _channelidx])
        # end
    end
    if channelids isa Number
        channelids = [channelids]
    end
    X = ToolsArray(lfp, (𝑡(timedata), Chan(channelids)))
    return X
end

# """
# Get the lfp data for a probe, providing *indices* for channels and times. See function below for indexing by channel ids and time values/intervals
# """
# function _getlfp(session::AbstractSession, probeid::Int; channelidxs=1:length(getlfpchannels(session, probeid)), timeidxs=1:length(getlfptimes(session, probeid)))
#     @assert(any(getprobeids(session) .== probeid), "Probe $probeid does not belong to session $(getid(session))")
#     @assert(subset(getprobes(session), :id=>ByRow(==(probeid)))[!, :has_lfp_data][1], @error "Probe $probeid does not have LFP data")
#     path = getlfppath(session, probeid)
#     if !isfile(path)
#         downloadlfp(session, probeid)
#     end
#     if !isfile(path)
#         @error "LFP data did not download correctly"
#     end

#     f = h5open(path)

#     timedata = getlfptimes(session, probeid, timeidxs)
#     dopermute = true
#     channelids = getlfpchannels(session, probeid)
#     channelids = channelids[channelidxs]
#     if (channelidxs isa Union{Int64, AbstractRange{Int64}}) & (timeidxs isa Union{Int64, AbstractRange{Int64}}) # Can use HDF5 slicing
#         lfp = f["acquisition"][splitext(basename(path))[1]][splitext(basename(path))[1]*"_data"]["data"][channelidxs, timeidxs]
#     elseif timeidxs isa Union{Int64, AbstractRange{Int64}}
#         lfp = [f["acquisition"][splitext(basename(path))[1]][splitext(basename(path))[1]*"_data"]["data"][i, timeidxs] for i ∈ channelidxs]
#         lfp = hcat(lfp...)
#         dopermute = false
#     else
#         lfp = read(f["acquisition"][splitext(basename(path))[1]][splitext(basename(path))[1]*"_data"]["data"])
#         lfp = lfp[channelidxs, timeidxs]
#     end
#     if lfp isa Vector
#        lfp = reshape(lfp, 1, length(lfp))
#     end
#     if channelids isa Number
#         channelids = [channelids]
#     end
#     if dopermute
#         lfp = permutedims(lfp, reverse(1:ndims(lfp)))
#     end
#     X = ToolsArray(lfp, (𝑡(timedata),  Chan(channelids)))
#     close(f)
#     return X
# end

using NWBStream
export AbstractNWBSession, NWBSession, S3Session

sessionfromnwb(S) = behavior_ecephys_session.BehaviorEcephysSession.from_nwb(S)

abstract type AbstractNWBSession <: AbstractSession end

struct NWBSession <: AbstractNWBSession
    file::Any
end
mutable struct S3Session <: AbstractNWBSession
    url::Any
    file::Any
    io::Any
    pyObject::Any
    function S3Session(url::String, file, io, pyObject = ())
        S = new(url, file, io, pyObject)
        f(S::S3Session) = @async s3close(S.io)
        finalizer(f, S)
    end
end
initialize!(S::AbstractNWBSession) = (S.pyObject = sessionfromnwb(S.file); nothing) # Can take absolutely forever since dataframes suck
S3Session(url::String) = (S = S3Session(url, s3open(url)...); initialize!(S); S)

getid(S::AbstractNWBSession) = pyconvert(Int, getfile(S).identifier |> string |> Meta.parse)
# getprobes(S::AbstractNWBSession) = Dict(getfile(S).ec_electrode_groups)
# getprobeids(S::AbstractNWBSession) = Dict(first(s) => pyconvert(Int, last(s).probe_id) for s in getprobes(S))
# getprobefiles(S::AbstractNWBSession) = Dict(getfile(S).ec_electrode_groups)
# getchannels() = manifest["project_metadata"]["channels.csv"] |> url2df

function getfile(S::AbstractNWBSession)
    return S.file
end

function getlfpchannels(session::S3Session)
    f = getfile(session)
    channels = f["general"]["extracellular_ephys"]["electrodes"]["id"][:]
    close(f)
    return channels
end

# getprobes(S::AbstractNWBSession) = S.pyObject.get_probes_obj().to_dataframe() |> PyPandasDataFrame |> DataFrame
function getsessionunits(session::AbstractNWBSession)
    units = session.pyObject.get_units() |> PyPandasDataFrame |> DataFrame
end

function py2df(p)
    df = p |> PyPandasDataFrame
    df = df |> DataFrame
    # if hasfield(p, :index)
    index = p.index
    index = pyconvert(Vector{Int64}, index)
    if string(p.index.name) == "None"
        df[!, "Index"] = index
    else
        df[!, string(p.index.name)] = index
    end
    # Manual conversions
    for s in names(df)
        if eltype(df[:, s]) <: Union{Missing, AbstractArray{Bool, 0}} # && all(length.(df[:, s]) .< 2)
            df[!, s] = [ismissing(i) ? missing : pyconvert(Bool, i[1]) for i in df[:, s]]
        end
        if eltype(df[:, s]) <: PythonCall.Wrap.PyArray
            df[!, s] = Array.(df[:, s])
        end
    end
    # end
    return df[!, [end, (1:(size(df, 2) - 1))...]]
end

# getid(S::AbstractNWBSession) = pyconvert(Int, S.pyObject.behavior_ecephys_session_id)
getprobes(S::AbstractNWBSession) = py2df(S.pyObject.probes)
getchannels(S::AbstractNWBSession) = py2df(S.pyObject.get_channels())

Base.Dict(p::Py) = pyconvert(Dict, p)
function Base.Dict(d::Dict{T, <:PythonCall.PyArray}) where {T}
    return Dict(k => pyconvert(Array{eltype(v)}, v) for (k, v) in d)
end

function getprobefile(session::AbstractNWBSession, name::AbstractString)
    files = getprobefiles(session)
    return files["probe_$(name)_lfp.nwb"]
end

function getprobefile(session::AbstractNWBSession, probeid::Int)
    getprobefile(session, first(subset(getprobes(session), :id => ByRow(==(probeid))).name))
end

function getlfpchannels(session::AbstractNWBSession, probeid::Int)
    f = getprobefile(session, probeid) |> NWBStream.s3open |> first
    _lfp = Dict(f.acquisition)["probe_1108501239_lfp_data"]
    channelids = pyconvert(Vector{Int64}, _lfp.electrodes.to_dataframe().index.values)
end

function getlfptimes(session::AbstractNWBSession, probeid::Int)
    f = getprobefile(session, probeid) |> NWBStream.s3open |> first
    _lfp = Dict(f.acquisition)["probe_1108501239_lfp_data"]
    timedata = pyconvert(Vector{Float32}, _lfp.timestamps)
end
function getlfptimes(session::AbstractNWBSession, probeid::Int, idxs)
    d = diff(idxs)
    if all(d[1] .== d)
        idxs = @py slice(idxs[0] - 1, idxs[-1] - 1, d[1])
    end
    f = getprobefile(session, probeid) |> NWBStream.s3open |> first
    _lfp = Dict(f.acquisition)["probe_1108501239_lfp_data"]
    timedata = _lfp.timestamps
    timedata = pyconvert(Vector{Float32}, timedata[idxs])
end
function getlfptimes(session::AbstractNWBSession, probeid::Int, i::Interval)
    f = getprobefile(session, probeid) |> NWBStream.s3open |> first
    _lfp = Dict(f.acquisition)["probe_1108501239_lfp_data"]
    timedata = _lfp.timestamps
    # infer the timestep
    dt = mean(diff(pyconvert(Vector, timedata[0:1000])))
    putativerange = pyconvert(Float64, timedata[0]):dt:pyconvert(Float64, timedata[-1] + dt)
    idxs = findall(putativerange .∈ (i,))
    idxs = vcat(idxs[1] .- (100:-1:1), idxs, idxs[end] .+ (1:100))
    idxs = idxs[idxs .> 0]
    idxs = idxs[idxs .< pyconvert(Int, timedata.len())]
    ts = getlfptimes(session::AbstractNWBSession, probeid::Int, idxs)
    ts = ts[ts .∈ (i,)]
end

function getepochs(S::AbstractNWBSession)
    df = S.pyObject.stimulus_presentations.groupby("stimulus_block").head(1) |> py2df
    df.stop_time = df.end_time
    return df
end

function _getlfp(session::AbstractNWBSession, probeid::Int;
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
    X = DimArray(lfp, (Ti(timedata), Dim{:channel}(channelids)))
    return X
end

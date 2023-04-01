using NWBS3
export AbstractNWBSession, NWBSession, S3Session

sessionfromnwb(S) = behavior_ecephys_session.BehaviorEcephysSession.from_nwb(S)

abstract type AbstractNWBSession <: AbstractSession end

struct NWBSession <: AbstractNWBSession
    file
end
mutable struct S3Session <: AbstractNWBSession
    url
    file
    io
    pyObject
    function S3Session(url::String, file, io, pyObject=())
        S = new(url, file, io, pyObject)
        f(S::S3Session) = @async s3close(S.io)
        finalizer(f, S)
    end
end
initialize!(S::AbstractNWBSession) = (S.pyObject = sessionfromnwb(S.file); nothing) # Can take absolutely forever since dataframes suck
S3Session(url::String) = S3Session(url, s3open(url)...)

getid(S::AbstractNWBSession) = getfile(S).identifier |> string |> Meta.parse
# getprobes(S::AbstractNWBSession) = Dict(getfile(S).ec_electrode_groups)
# getprobeids(S::AbstractNWBSession) = Dict(first(s) => pyconvert(Int, last(s).probe_id) for s in getprobes(S))
# getprobefiles(S::AbstractNWBSession) = Dict(getfile(S).ec_electrode_groups)
# getchannels() = manifest["project_metadata"]["channels.csv"] |> url2df

function getfile(S::AbstractNWBSession)
    return S.file
end

function getlfpchannels(session::S3Session, probeid)
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
        df[!, string(p.index.name)] = index
    # end
    return df[!, [end, (1:size(df, 2)-1)...]]
end

getid(S::AbstractNWBSession) = S.pyObject.behavior_session_id
getprobes(S::AbstractNWBSession) = py2df(S.pyObject.probes)
getchannels(S::AbstractNWBSession) = py2df(S.pyObject.get_channels())
getepochs(S::AbstractNWBSession) = S.pyObject.stimulus_presentations.groupby("stimulus_block").head(1) |> py2df

Base.Dict(p::Py) = pyconvert(Dict, p)

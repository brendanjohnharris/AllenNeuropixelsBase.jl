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
    function S3Session(url::String, file, io, pyobject=())
        S = new(url, file, io, pyobject)
        f(S::S3Session) = @async s3close(S.io)
        finalizer(f, S)
    end
end
initialize!(S::AbstractNWBSession) = (S.pyobject = sessionfromnwb(S.file); nothing) # Can take absolutely forever since dataframes suck
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


Base.Dict(p::Py) = pyconvert(Dict, p)

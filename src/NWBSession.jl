using NWBS3
export AbstractNWBSession, NWBSession, S3Session

abstract type AbstractNWBSession <: AbstractSession end

struct NWBSession <: AbstractNWBSession
    file
end
struct S3Session <: AbstractNWBSession
    url
    file
    S3Session(url::String, file=s3open(url, "rb")) = New(url, file)
end

function getfile(S::AbstractNWBSession)
    return S.file
end



function getlfpchannels(session::S3Session, probeid)
    f = getfile(session)
    channels = f["general"]["extracellular_ephys"]["electrodes"]["id"][:]
    close(f)
    return channels
end

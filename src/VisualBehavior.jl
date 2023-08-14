module VisualBehavior
# using ..AllenNeuropixelsBase
import ..AllenNeuropixelsBase as ANB
using PythonCall
using NWBStream
using DataFrames
using Downloads
using JSON
using HDF5
import NWBStream.url2df
import ..AllenNeuropixelsBase: AbstractNWBSession, AbstractSession, py2df

const visualbehavior = "https://visual-behavior-neuropixels-data.s3.us-west-2.amazonaws.com/"
const visualbehaviorbehavior = visualbehavior * "visual-behavior-neuropixels/"

"""
    getmanifesturl(manifest_version="0.4.0")

Return the URL for the neuropixels visual behavior project manifest file hosted on a given `hostname` server and version `manifest_version`.

## Arguments
- `manifest_version::String`: The version of the manifest file. Default value is `"0.4.0"`.

Other manifest version numbers can be identified here: https://s3.console.aws.amazon.com/s3/buckets/visual-behavior-neuropixels-data?prefix=visual-behavior-neuropixels%2Fmanifests%2F&region=us-west-2

"""
function getmanifesturl(manifest_version="0.4.0")
    object_key = "manifests/visual-behavior-neuropixels_project_manifest_v$manifest_version.json"
    return visualbehaviorbehavior * object_key
end

function getmanifest(args...)
    hostname = getmanifesturl(args...)
    data = String(take!(Downloads.download(hostname, IOBuffer())))
    JSON.parse(data)
end

"""
    filetree(paths::Vector{String})

Creates a dictionary-like object representing the directory and file structure from a list of file paths.

## Arguments:
 * paths: A vector of file paths to include in the tree.

## Returns:
 * dict_tree : A dictionary-like object representing the directory and file structure.
"""
function filetree(paths::Vector{String})
    dict_tree = Dict{String,Any}()
    for path in paths
        current_dict = dict_tree
        components = splitpath(path)
        for component in components[1:end-1]
            if !haskey(current_dict, component)
                current_dict[component] = Dict{String,Any}()
            end
            current_dict = current_dict[component]
        end
        current_dict[components[end]] = path
    end
    return dict_tree
end


"""
     datapaths()

Get the paths to the visual behavior neuropixels data on Amazon S3 as a tree-like dictionary.
"""
function datapaths()
    data = getmanifest()
    ks = merge([data[k] for k in keys(data) if data[k] isa Dict]...)
    ks = getindex.([ks[k] isa Dict ? ks[k] : Dict(k => ks[k]) for k in keys(ks)], ["url"])
    tree = filetree(String.(ks))
    return tree["https:"]["visual-behavior-neuropixels-data.s3.us-west-2.amazonaws.com"]["visual-behavior-neuropixels"]
end

const manifest = datapaths()


struct Session <: AbstractSession
    pyObject
    behavior_pyObject
    function Session(pyObject)
        behavior_pyObject = ANB.behaviorcache().get_behavior_session(pyObject.behavior_session_id)
        new(pyObject, behavior_pyObject)
    end
end
function Session(session_id::Int; kwargs...)
    ANB.getsessiondata(session_id; kwargs...) |> Session
end
Session(; params...) = Session(params[:sessionid]);


getprobes() = manifest["project_metadata"]["probes.csv"] |> url2df
getchannels() = manifest["project_metadata"]["channels.csv"] |> url2df
ANB.getid(S::Session) = pyconvert(Int, S.pyObject.metadata["ecephys_session_id"])
getbehaviorid(S::Session) = pyconvert(Int, S.pyObject.behavior_session_id)

function ANB.getchannels(S::Session, args...; kwargs...)
    df = py2df(S.pyObject.get_channels())
    ANB.getchannels(df, args...; kwargs...)
end

const sources = Dict(
    :Allen_neuropixels_visual_behavior => "https://visual-behavior-neuropixels-data.s3.us-west-2.amazonaws.com"
)

function ANB.getlfppath(session::Session, probeid)
    probe = ANB.getprobeobj(session, probeid)
    @assert pyconvert(Int, probe.id) == probeid
    name = probe.name
    path = pyconvert(String, probe._lfp_meta.lfp_csd_filepath()._str)
end

function ANB.isinvalidtime(session::Session, probeids=getprobeids(session), times=NaN)
    # paths = ANB.getlfppath.((session,), probeids)
    # for path in paths
    #     f = h5open(path)
    # end
    false # no invalid times for VBN as far as I can see
end

function ANB.getlfptimes(session::Session, probeid)
    path = ANB.getlfppath(session, probeid)
    if !isfile(path)
        @error "Probe lfp data file does not exist. This may be becuase you have not downloaded the probe data. Try using `getlfp`"
    end

    f = h5open(path)
    stem = ANB.resolvenwblfp(session, probeid)
    times = f["acquisition"][stem][stem*"_data"]["timestamps"][:]
    close(f)
    return times
end

S3Session(session_id::Int, args...; kwargs...) = ANB.S3Session(getsessionfile(session_id), args...; kwargs...)

"""
    getsessionfiles()

Return a dictionary mapping session IDs to their corresponding file paths.
"""
function getsessionfiles()
    files = manifest["behavior_ecephys_sessions"]
    sessionids = keys(files)
    D = Dict()
    for s in sessionids
        subfiles = files[s]
        subfile = findfirst(subfiles) do f
            occursin("ecephys_session_$s", f)
        end
        push!(D, s => subfiles[subfile])
    end
    return D
end

getsessionfile(session_id::Int, args...) = getsessionfiles(args...)["$session_id"]

function getsessiontable()
    st = manifest["project_metadata"]["ecephys_sessions.csv"] |> url2df
    st.structure_acronyms = [replace.(split(a, ','), r"[\[\]\'\s]"=>"") for a in st.structure_acronyms]
    st.ecephys_structure_acronyms = st.structure_acronyms
    return st
end
getunits() = manifest["project_metadata"]["units.csv"] |> url2df
function getunitanalysismetricsbysessiontype()
    session_table = ANB.VisualBehavior.getsessiontable()
    units = getunits()
    sessionmetrics = innerjoin(session_table, metrics, on=:ecephys_session_id)
end

function getprobes(session::AbstractNWBSession)
    probes = VisualBehavior.getprobes()
    return subset(probes, :ecephys_session_id => ByRow(==(ANB.getid(session))))
end

function getprobefiles()
    files = manifest["behavior_ecephys_sessions"]
    sessionids = keys(files)
    D = deepcopy(files)
    for s in sessionids
        subfiles = D[s]
        subfile = findfirst(subfiles) do f
            occursin("ecephys_session_$s", f)
        end
        delete!(subfiles, subfile)
    end
    return D
end
getprobefiles(S::AbstractNWBSession) = getprobefiles()[string(ANB.getid(S))]

function getprobefile(session::AbstractNWBSession, name::AbstractString)
    files = getprobefiles(session)
    return files["probe_$(name)_lfp.nwb"]
end

getprobefile(session::AbstractNWBSession, probeid::Int) = getprobefile(session, first(subset(getprobes(session), :ecephys_probe_id => ByRow(==(probeid))).name))

function ANB.getepochs(S::Session)
    df = S.pyObject.stimulus_presentations.groupby("stimulus_block").head(1) |> py2df
    df.stop_time = df.end_time
    return df
end

# * The goal is now to copy all of the convenience functions in AllenAttention.jl to here...

# Will just need to overload getlfp and getspikes, since we can use the allensdk to get all the session metadata and such.
# Will need VisualBehaviorNeuropixelsProjectCache
# from allensdk.brain_observatory.behavior.behavior_project_cache import VisualBehaviorNeuropixelsProjectCache
# See from allensdk.brain_observatory.behavior.behavior_project_cache import VisualBehaviorNeuropixelsProjectCache
# All the api above should go into allen neuropixels. This here is a small package.

# Session(sessionid::Int) = sessionid |> getsessionfile |> ANB.Session

# function getlfp(session::S3Session, probeid::Int; channels=getlfpchannels(session, probeid), times=ClosedInterval(extrema(getlfptimes(session, probeid))...), inbrain=false)
# end

# VisualBehaviorNeuropixelsProjectCache = ANB.behavior_project_cache.behavior_neuropixels_project_cache.VisualBehaviorNeuropixelsProjectCache
# VisualBehaviorNeuropixelsProjectCache.from_s3_cache(mktempdir())

function ANB.getchannels(S::Session)
    s = py2df(S.pyObject.get_channels())
end


function getunitanalysismetrics()
    session_table = ANB.VisualBehavior.getsessiontable()
    metrics = ANB.behaviorcache().get_unit_table() |> py2df
    metrics = outerjoin(metrics, session_table; on=:ecephys_session_id)
end

function ANB.getunitmetrics(session::Session)
    str = session.pyObject.get_units()
    py2df(str)
end

function ANB.getprobeobj(S::Session, probeid)
    probes = S.pyObject.get_probes_obj()
    probeids = [pyconvert(Int, probe.id) for probe in probes.probes]
    thisprobe = findfirst(probeids .== probeid)-1
    probe = probes.probes[thisprobe]
end


end
export VisualBehavior
getprobefiles(S::AbstractNWBSession; dataset=VisualBehavior) = dataset.getprobefiles(S)

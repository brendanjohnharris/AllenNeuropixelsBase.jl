module VisualBehavior
using ..AllenNeuropixelsBase
import ..AllenNeuropixelsBase as ANB
using NWBS3
using DataFrames
using Downloads
using JSON
import NWBS3.url2df

const visualbehavior = "https://visual-behavior-neuropixels-data.s3.us-west-2.amazonaws.com/"
const visualbehaviorbehaviour = visualbehavior * "visual-behavior-neuropixels/"

"""
    getmanifesturl(manifest_version="0.4.0")

Return the URL for the neuropixels visual behavior project manifest file hosted on a given `hostname` server and version `manifest_version`.

## Arguments
- `manifest_version::String`: The version of the manifest file. Default value is `"0.4.0"`.

Other manifest version numbers can be identified here: https://s3.console.aws.amazon.com/s3/buckets/visual-behavior-neuropixels-data?prefix=visual-behavior-neuropixels%2Fmanifests%2F&region=us-west-2

"""
function getmanifesturl(manifest_version="0.4.0")
    object_key = "manifests/visual-behavior-neuropixels_project_manifest_v$manifest_version.json"
    return visualbehaviorbehaviour * object_key
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

getprobes() = manifest["project_metadata"]["probes.csv"] |> url2df
getchannels() = manifest["project_metadata"]["channels.csv"] |> url2df

const sources = Dict(
    :Allen_neuropixels_visual_behaviour => "https://visual-behavior-neuropixels-data.s3.us-west-2.amazonaws.com"
)



S3Session(session_id::Int, args...; kwargs...) = S3Session(getsessionfile(session_id), args...; kwargs...)


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

getsessiontable() = manifest["project_metadata"]["ecephys_sessions.csv"] |> url2df


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



# * The goal is now to copy all of the convenience functions in AllenAttention.jl to here...

# Will just need to overload getlfp and getspikes, since we can use the allensdk to get all the session metadata and such.
# Will need VisualBehaviorNeuropixelsProjectCache
# from allensdk.brain_observatory.behavior.behavior_project_cache import VisualBehaviorNeuropixelsProjectCache
# See from allensdk.brain_observatory.behavior.behavior_project_cache import VisualBehaviorNeuropixelsProjectCache
# All the api above should go into allen neuropixels. This here is a small package.

Session(sessionid::Int) = sessionid |> getsessionfile |> ANB.S3Session

# function getlfp(session::S3Session, probeid::Int; channels=getlfpchannels(session, probeid), times=ClosedInterval(extrema(getlfptimes(session, probeid))...), inbrain=false)
# end

# VisualBehaviorNeuropixelsProjectCache = ANB.behavior_project_cache.behavior_neuropixels_project_cache.VisualBehaviorNeuropixelsProjectCache
# VisualBehaviorNeuropixelsProjectCache.from_s3_cache(mktempdir())


end
export VisualBehavior

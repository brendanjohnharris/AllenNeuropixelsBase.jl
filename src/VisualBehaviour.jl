module VisualBehaviour
import ..AllenNeuropixelsBase as ANB
using NWBS3
using DataFrames
using Downloads
using JSON
import NWBS3.url2df

const visualbehaviour = "https://visual-behavior-neuropixels-data.s3.us-west-2.amazonaws.com/"
const visualbehaviourbehaviour = visualbehaviour * "visual-behavior-neuropixels/"


function getmanifesturl(manifest_version="0.4.0")
    object_key = "manifests/visual-behavior-neuropixels_project_manifest_v$manifest_version.json"
    return visualbehaviourbehaviour * object_key
end

function getmanifest(args...)
    hostname = getmanifesturl(args...)
    data = String(take!(Downloads.download(hostname, IOBuffer())))
    JSON.parse(data)
end

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

function datapaths()
    data = getmanifest()
    ks = merge([data[k] for k in keys(data) if data[k] isa Dict]...)
    ks = getindex.([ks[k] isa Dict ? ks[k] : Dict(k => ks[k]) for k in keys(ks)], ["url"])
    tree = filetree(String.(ks))
    return tree["https:"]["visual-behavior-neuropixels-data.s3.us-west-2.amazonaws.com"]["visual-behavior-neuropixels"]
end

const manifest = datapaths()

const sources = Dict(
    :Allen_neuropixels_visual_behaviour => "https://visual-behavior-neuropixels-data.s3.us-west-2.amazonaws.com"
)

"""
getmanifesturl(hostname=visualbehaviour, manifest_version="0.1.0")

Return the URL for the neuropixels visual behavior project manifest file hosted on a given `hostname` server and version `manifest_version`.

## Arguments
- `hostname::String`: The hostname of the server where the manifest file is hosted. Default value is `visualbehaviour`.
- `manifest_version::String`: The version of the manifest file. Default value is `"0.4.0"`.

Other manifest version numbers can be identified here: https://s3.console.aws.amazon.com/s3/buckets/visual-behavior-neuropixels-data?prefix=visual-behavior-neuropixels%2Fmanifests%2F&region=us-west-2

"""
function getmanifesturl(hostname=visualbehaviour, manifest_version="0.4.0")
    object_key = "visual-behavior-neuropixels/manifests/visual-behavior-neuropixels_project_manifest_v$manifest_version.json"
    return hostname * object_key
end



S3Session(session_id::Int, args...; kwargs...) = S3Session(getsessionfile(session_id), args...; kwargs...)

getprobes() = manifest["project_metadata"]["probes.csv"] |> url2df

getchannels() = manifest["project_metadata"]["channels.csv"] |> url2df

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
            occursin(s, f)
        end
        push!(D, s => subfiles[subfile])
    end
    return D
end

getsessionfile(session_id::Int, args...) = getsessionfiles(args...)["ecephys_session_$session_id.nwb"]

getsessiontable() = manifest["project_metadata"]["ecephys_sessions.csv"] |> url2df


# * The goal is now to copy all of the convenience functions in AllenAttention.jl to here...

# Will just need to overload getlfp and getspikes, since we can use the allensdk to get all the session metadata and such.
# Will need VisualBehaviorNeuropixelsProjectCache
# from allensdk.brain_observatory.behavior.behavior_project_cache import VisualBehaviorNeuropixelsProjectCache
# See from allensdk.brain_observatory.behavior.behavior_project_cache import VisualBehaviorNeuropixelsProjectCache
# All the api above should go into allen neuropixels. This here is a small package.

Session(sessionid::Int) = sessionid |> getsessionfile |> ANB.S3Session

# function getlfp(session::S3Session, probeid::Int; channels=getlfpchannels(session, probeid), times=ClosedInterval(extrema(getlfptimes(session, probeid))...), inbrain=false)
# end

end

module VisualBehaviour
using NWBS3
using DataFrames
using Downloads
using JSON
import NWBS3.url2df

const visualbehaviour = "https://visual-behavior-neuropixels-data.s3.us-west-2.amazonaws.com/"
const visualbehaviourbehaviour = visualbehaviour * "visual-behavior-neuropixels/"


function getmanifesturl(manifest_version="0.1.0")
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


function getmanifesturl(hostname=visualbehaviour, manifest_version="0.1.0")
    object_key = "visual-behavior-neuropixels/manifests/visual-behavior-neuropixels_project_manifest_v$manifest_version.json"
    return hostname * object_key
end



NWBS3.S3Session(session_id::Int, args...; kwargs...) = NWBS3.S3Session(getsessionfile(session_id), args...; kwargs...)

getprobes() = manifest["project_metadata"]["probes.csv"] |> url2df

getchannels() = manifest["project_metadata"]["channels.csv"] |> url2df


getsessionfiles() = manifest["ecephys_sessions"]

getsessionfile(session_id::Int, args...) = getsessionfiles(args...)["ecephys_session_$session_id.nwb"]

getsessiontable() = manifest["project_metadata"]["ecephys_sessions.csv"] |> url2df


# * The goal is now to copy all of the convenience functions in AllenAttention.jl to here...

# Will just need to overload getlfp and getspikes, since we can use the allensdk to get all the session metadata and such.
# Will need VisualBehaviorNeuropixelsProjectCache
# from allensdk.brain_observatory.behavior.behavior_project_cache import VisualBehaviorNeuropixelsProjectCache
# See from allensdk.brain_observatory.behavior.behavior_project_cache import VisualBehaviorNeuropixelsProjectCache
# All the api above should go into allen neuropixels. This here is a small package.
end

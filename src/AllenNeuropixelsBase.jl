module AllenNeuropixelsBase
using PythonCall
using DataFrames
using CSV
using TimeseriesTools
using Preferences
using AllenSDK
using Downloads

const brain_observatory = PythonCall.pynew()
const stimulusmapping = PythonCall.pynew()
const ecephys_project_cache = PythonCall.pynew()
const ecephys_project_api = PythonCall.pynew()
const ephys_features = PythonCall.pynew()
const brain_observatory_cache = PythonCall.pynew()
const stimulus_info = PythonCall.pynew()
const mouse_connectivity_cache = PythonCall.pynew()
const ontologies_api = PythonCall.pynew()
const reference_space_cache = PythonCall.pynew()
const reference_space = PythonCall.pynew()
const nwb_api = PythonCall.pynew()
const behavior_project_cache = PythonCall.pynew()
const behavior_ecephys_session = PythonCall.pynew()
const _behaviorcache = PythonCall.pynew()
export allensdk, brain_observatory, ecephys, ecephys_project_cache,
       mouse_connectivity_cache, ontologies_api, reference_space_cache, reference_space,
       behavior_ecephys_session, behavior_project_cache

function setdatadir(datadir::String)
    @set_preferences!("datadir"=>datadir)
    @info("New default datadir set; restart your Julia session for this change to take effect")
end
const datadir = replace(@load_preference("datadir",
                                         joinpath(pkgdir(AllenNeuropixelsBase), "data/")),
                        "\\" => "/")
const ecephysmanifest = replace(joinpath(datadir, "Ecephys", "manifest.json"), "\\" => "/")
const behaviormanifest = replace(joinpath(datadir, "Behavior"),
                                 "\\" => "/")
const brainobservatorymanifest = replace(joinpath(datadir, "BrainObservatory",
                                                  "manifest.json"), "\\" => "/")
const mouseconnectivitymanifest = replace(joinpath(datadir, "MouseConnectivity",
                                                   "manifest.json"), "\\" => "/")
const referencespacemanifest = replace(joinpath(datadir, "ReferenceSpace", "manifest.json"),
                                       "\\" => "/")
export setdatadir, datadir, ecephysmanifest, brainobservatorymanifest,
       mouseconnectivitymanifest, referencespacemanifest, behaviormanifest

const streamlinepath = abspath(referencespacemanifest, "../laplacian_10.nrrd")

function __init__()
    PythonCall.pycopy!(brain_observatory, pyimport("allensdk.brain_observatory"))
    PythonCall.pycopy!(stimulus_info, pyimport("allensdk.brain_observatory.stimulus_info"))
    PythonCall.pycopy!(stimulusmapping,
                       pyimport("allensdk.brain_observatory.ecephys.stimulus_analysis.receptive_field_mapping"))
    PythonCall.pycopy!(ecephys_project_cache,
                       pyimport("allensdk.brain_observatory.ecephys.ecephys_project_cache"))
    PythonCall.pycopy!(ecephys_project_api,
                       pyimport("allensdk.brain_observatory.ecephys.ecephys_project_api"))
    PythonCall.pycopy!(ephys_features, pyimport("allensdk.ephys.ephys_features"))
    PythonCall.pycopy!(brain_observatory_cache,
                       pyimport("allensdk.core.brain_observatory_cache"))
    PythonCall.pycopy!(mouse_connectivity_cache,
                       pyimport("allensdk.core.mouse_connectivity_cache"))
    PythonCall.pycopy!(ontologies_api, pyimport("allensdk.api.queries.ontologies_api"))
    PythonCall.pycopy!(reference_space_cache,
                       pyimport("allensdk.core.reference_space_cache"))
    PythonCall.pycopy!(reference_space, pyimport("allensdk.core.reference_space"))
    PythonCall.pycopy!(nwb_api, pyimport("allensdk.brain_observatory.nwb.nwb_api"))
    PythonCall.pycopy!(behavior_project_cache,
                       pyimport("allensdk.brain_observatory.behavior.behavior_project_cache"))
    PythonCall.pycopy!(behavior_ecephys_session,
                       pyimport("allensdk.brain_observatory.ecephys.behavior_ecephys_session"))

    ecephys_project_cache.EcephysProjectCache.from_warehouse(manifest = ecephysmanifest)

    if !isfile(streamlinepath)
        mkpath(dirname(streamlinepath))
        @info "Downloading streamline data to $streamlinepath, this may take a few minutes"
        Downloads.download("https://www.dropbox.com/sh/7me5sdmyt5wcxwu/AACFY9PQ6c79AiTsP8naYZUoa/laplacian_10.nrrd?dl=1",
                           streamlinepath)
        @info "Download streamline data with status $(isfile(streamlinepath))"
    end

    PythonCall.pycopy!(_behaviorcache, __behaviorcache())
end

"""
    loaddataframe(file, dir=datadir)::DataFrame

Load a CSV file into a DataFrame.

Arguments:

- `file`: A string representing the name of the CSV file to be loaded.
- `dir` (optional): A string representing the directory containing the CSV file. Default is `datadir`.

Returns:

A `DataFrame` object containing the contents of the CSV file.

Example:
```julia
using DataFrames
df = loaddataframe("mydata.csv", "/path/to/my/data/")
```
"""
function loaddataframe(file, dir = datadir)::DataFrame
    CSV.File(abspath(dir, file)) |> DataFrame
end
export convertdataframe

include("./EcephysCache.jl")
include("./BrainObservatory.jl")
# include("./HybridSession.jl")
include("./SparseDimArray.jl")
include("./MouseConnectivityCache.jl")
include("./Ontologies.jl")
include("./ReferenceSpace.jl")
include("./NWBSession.jl")
include("./LFP.jl")
include("./SpikeBand.jl")
include("./Behaviour.jl")
include("./VisualBehavior.jl")

end

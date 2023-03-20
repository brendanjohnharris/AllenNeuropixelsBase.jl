module AllenNeuropixelsBase
using PythonCall
using DataFrames
using CSV
using Preferences

const pynwb = PythonCall.pynew()
const allensdk = PythonCall.pynew()
const brain_observatory = PythonCall.pynew()
const ecephys = PythonCall.pynew()
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
export allensdk, brain_observatory, ecephys, ecephys_project_cache, mouse_connectivity_cache, ontologies_api, reference_space_cache, reference_space


function __init__()
    PythonCall.pycopy!(pynwb, pyimport("pynwb"))
    PythonCall.pycopy!(allensdk, pyimport("allensdk"))
    PythonCall.pycopy!(brain_observatory, pyimport("allensdk.brain_observatory"))
#     PythonCall.pycopy!(stimulus_info, pyimport("allensdk.brain_observatory.stimulus_info"))
#     PythonCall.pycopy!(ecephys, pyimport("allensdk.brain_observatory.ecephys"))
#     PythonCall.pycopy!(stimulusmapping, pyimport("allensdk.brain_observatory.ecephys.stimulus_analysis.receptive_field_mapping"))
#     PythonCall.pycopy!(ecephys_project_cache, pyimport("allensdk.brain_observatory.ecephys.ecephys_project_cache"))
#     PythonCall.pycopy!(ecephys_project_api, pyimport("allensdk.brain_observatory.ecephys.ecephys_project_api"))
#     PythonCall.pycopy!(ephys_features, pyimport("allensdk.ephys.ephys_features"))
#     PythonCall.pycopy!(brain_observatory_cache, pyimport("allensdk.core.brain_observatory_cache"))
#     PythonCall.pycopy!(mouse_connectivity_cache, pyimport("allensdk.core.mouse_connectivity_cache"))
#     PythonCall.pycopy!(ontologies_api, pyimport("allensdk.api.queries.ontologies_api"))
#     PythonCall.pycopy!(reference_space_cache, pyimport("allensdk.core.reference_space_cache"))
#     PythonCall.pycopy!(reference_space, pyimport("allensdk.core.reference_space"))

#     ecephys_project_cache.EcephysProjectCache.from_warehouse(manifest=ecephysmanifest)
end

function setdatadir(datadir::String)
    @set_preferences!("datadir" => datadir)
    @info("New default datadir set; restart your Julia session for this change to take effect")
end
const datadir = replace(@load_preference("datadir", joinpath(pkgdir(AllenNeuropixelsBase), "data/")), "\\"=>"/")
const ecephysmanifest = replace(joinpath(datadir, "Ecephys", "manifest.json"), "\\"=>"/")
const brainobservatorymanifest = replace(joinpath(datadir, "BrainObservatory", "manifest.json"), "\\"=>"/")
const mouseconnectivitymanifest = replace(joinpath(datadir, "MouseConnectivity", "manifest.json"), "\\"=>"/")
const referencespacemanifest = replace(joinpath(datadir, "ReferenceSpace", "manifest.json"), "\\"=>"/")
export setdatadir, datadir, ecephysmanifest, brainobservatorymanifest, mouseconnectivitymanifest, referencespacemanifest


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
function loaddataframe(file, dir=datadir)::DataFrame
    CSV.File(abspath(dir, file)) |> DataFrame
end
export convertdataframe

# include("./EcephysCache.jl")
# include("./BrainObservatory.jl")
# include("./VisualBehaviour.jl")
# include("./SparseDimArray.jl")
# include("./LFP.jl")
# include("./SpikeBand.jl")
# include("./MouseConnectivityCache.jl")
# include("./Ontologies.jl")
# include("./ReferenceSpace.jl")
# include("./NWBSession.jl")
# include("./Behaviour.jl")

end

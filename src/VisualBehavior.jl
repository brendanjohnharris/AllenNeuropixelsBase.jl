function getchangetrials end
function getchangetrialtimeseries end
export getchangetrials, getchangetrialtimeseries, stimulustimeseries

module VisualBehavior
# using ..AllenNeuropixelsBase
import ..AllenNeuropixelsBase as ANB
using PythonCall
using NWBStream
using DataFrames
using Downloads
using JSON
using HDF5
using IntervalSets
using TimeseriesTools
import NWBStream.url2df
import ..AllenNeuropixelsBase: AbstractNWBSession, AbstractSession, py2df
import DataFrames.groupby

const visualbehavior = "https://visual-behavior-neuropixels-data.s3.us-west-2.amazonaws.com/"
const visualbehaviorbehavior = visualbehavior * "visual-behavior-neuropixels/"

"""
    getmanifesturl(manifest_version="0.5.0")

Return the URL for the neuropixels visual behavior project manifest file hosted on a given `hostname` server and version `manifest_version`.

## Arguments
- `manifest_version::String`: The version of the manifest file. Default value is `"0.5.0"`.

Other manifest version numbers can be identified here: https://s3.console.aws.amazon.com/s3/buckets/visual-behavior-neuropixels-data?prefix=visual-behavior-neuropixels%2Fmanifests%2F&region=us-west-2

"""
function getmanifesturl(manifest_version = "0.5.0")
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
    dict_tree = Dict{String, Any}()
    for path in paths
        current_dict = dict_tree
        components = splitpath(path)
        for component in components[1:(end - 1)]
            if !haskey(current_dict, component)
                current_dict[component] = Dict{String, Any}()
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
    pyObject::Any
    behavior_pyObject::Any
    function Session(pyObject)
        behavior_pyObject = ANB.behaviorcache().get_behavior_session(pyObject.behavior_session_id)
        new(pyObject, behavior_pyObject)
    end
end
function Session(session_id::Int; kwargs...)
    @debug "Constructing a session from id $session_id"
    ANB.getsessiondata(session_id; kwargs...) |> Session
end
Session(; params...) = Session(params[:sessionid])

getprobes() = manifest["project_metadata"]["probes.csv"] |> url2df
getchannels() = manifest["project_metadata"]["channels.csv"] |> url2df
ANB.getid(S::Session) = pyconvert(Int, S.pyObject.metadata["ecephys_session_id"])
getbehaviorid(S::Session) = pyconvert(Int, S.pyObject.behavior_session_id)

function ANB.getchannels(S::Session, args...; kwargs...)
    df = py2df(S.pyObject.get_channels())
    ANB.getchannels(df, args...; kwargs...)
end

const sources = Dict(:Allen_neuropixels_visual_behavior => "https://visual-behavior-neuropixels-data.s3.us-west-2.amazonaws.com")

function ANB.getlfppath(session::Session, probeid)
    probe = ANB.getprobeobj(session, probeid)
    @assert pyconvert(Int, probe.id) == probeid
    name = probe.name
    try
        path = pyconvert(String, probe._lfp_meta.lfp_csd_filepath()._str)
    catch # Probe has no file name for some reason. Make a hacky guess
        @warn "Malformed probe object. Guessing path"
        path = joinpath(ANB.behaviormanifest, "visual-behavior-neuropixels-0.5.0",
                        "behavior_ecephys_sessions",
                        string(pyconvert(Int,
                                         session.pyObject.metadata["ecephys_session_id"])),
                        "probe_$(name)_lfp.nwb")
    end
end

function ANB.isinvalidtime(session::Session, probeids = getprobeids(session), times = NaN)
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
    times = f["acquisition"][stem][stem * "_data"]["timestamps"][:]
    close(f)
    return times
end

function ANB.getperformancemetrics(session::Session)
    pyconvert(Dict, session.behavior_pyObject.get_performance_metrics())
end

function S3Session(session_id::Int, args...; kwargs...)
    ANB.S3Session(getsessionfile(session_id), args...; kwargs...)
end

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
    st.structure_acronyms = [replace.(split(a, ','), r"[\[\]\'\s]" => "")
                             for a in st.structure_acronyms]
    st.ecephys_structure_acronyms = st.structure_acronyms
    return st
end

function getbehaviorsessiontable()
    st = manifest["project_metadata"]["behavior_sessions.csv"] |> url2df
    return st
end

getunits() = manifest["project_metadata"]["units.csv"] |> url2df
function getunitanalysismetricsbysessiontype()
    session_table = ANB.VisualBehavior.getsessiontable()
    units = getunits()
    sessionmetrics = innerjoin(session_table, metrics, on = :ecephys_session_id)
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

function getprobefile(session::AbstractNWBSession, probeid::Int)
    getprobefile(session,
                 first(subset(getprobes(session),
                              :ecephys_probe_id => ByRow(==(probeid))).name))
end

function ANB.getepochs(S::Session)
    df = S.pyObject.stimulus_presentations |> py2df
    df.interval = [a .. b for (a, b) in zip(df.start_time, df.end_time)]
    df = groupby(df, :stimulus_block)
    _intersect(x) = minimum(minimum.(x)) .. maximum(maximum.(x))
    df = DataFrames.combine(df, :interval => _intersect => :interval,
                            :active => unique => :active,
                            :stimulus_block => unique => :stimulus_block,
                            :stimulus_name => unique => :stimulus_name)
    df.start_time = minimum.(df.interval)
    df.end_time = maximum.(df.interval)
    df.duration = df.end_time - df.start_time
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

function getunitanalysismetrics(; filter_by_validity = nothing) # Extra kwarg for Visual Coding compat
    session_table = ANB.VisualBehavior.getsessiontable()
    metrics = ANB.behaviorcache().get_unit_table() |> py2df
    metrics = outerjoin(metrics, session_table; on = :ecephys_session_id)
end

function ANB.getunitanalysismetrics(session::Session; kwargs...)
    ANB.getunitmetrics(session; kwargs...)
end

function ANB.getunitmetrics(session::Session; filter_by_validity = false)
    if filter_by_validity
        @warn "Filtering by validity is not implemented for visual behavior neuropixels data"
    end
    str = session.pyObject.get_units()
    df = py2df(str)
    df.ecephys_unit_id = df.id
    df.ecephys_channel_id = df.peak_channel_id
    df
end

function ANB.getprobeobj(S::Session, probeid)
    probes = S.pyObject.get_probes_obj()
    probeids = [pyconvert(Int, probe.id) for probe in probes.probes]
    i = findfirst(probeids .== probeid)
    if isempty(i)
        @error "Probe id $probeid not found in session $(getid(S))"
    end
    thisprobe = i - 1
    probe = probes.probes[thisprobe]
end

## Behavior metrics
import AllenNeuropixelsBase: gettrials, getlicks, getrewards, geteyetracking,
                             getrunningspeed, getbehavior, getchangetrials,
                             getchangetrialtimeseries
function gettrials(S::Session)
    df = S.pyObject.trials |> py2df
    df.interval = [a .. b for (a, b) in zip(df.start_time, df.stop_time)]
    return df
end
getlicks(S::Session) = S.pyObject.licks |> py2df
getrewards(S::Session) = S.pyObject.rewards |> py2df
geteyetracking(S::Session) = S.pyObject.eye_tracking |> py2df
_getrunningspeed(S::Session) = S.pyObject.running_speed |> py2df

function _getbehavior(S::Session)
    fs = (gettrials, getlicks, getrewards, geteyetracking, getrunningspeed)
    dfs = [f(S) for f in fs]
end

function stimulustimeseries(session; blanks = true)
    df = ANB.getstimuli(session)
    x = TimeseriesTools.TimeSeries(df.start_time, df.image_name)
    y = TimeseriesTools.TimeSeries(df.end_time, df.image_name)
    if blanks
        xx = TimeseriesTools.TimeSeries(times(x) .- 1 / 1250, fill("blank", length(x)))
        yy = TimeseriesTools.TimeSeries(times(y) .+ 1 / 1250, fill("blank", length(y)))
        x = TimeseriesTools.interlace(x, xx)
        y = TimeseriesTools.interlace(y, yy)
    end
    TimeseriesTools.interlace(x, y)
end

function getchangetrials(session)
    df = gettrials(session)
    df = df[df.is_change, :]

    # Get corrected change time by matching stimulus frame with start time in stimulus table
    stimuli = ANB.getstimuli(session)
    stimuli.change_frame = stimuli.start_frame
    stimuli = innerjoin(df, stimuli, on = :change_frame, makeunique = true)
    df.change_time_with_display_delay = stimuli.start_time_1

    # Get response (lick) latency
    df.lick_latency = df.response_time - df.change_time_with_display_delay # ? Do we want to use the precise lick times from the licks dataframe?

    return df
end

function getchangetrialtimeseries(session)
    # Make an irregular time series of when the change trials occur
    trials = getchangetrials(session)

    t = trials.change_time_with_display_delay
    vars = [:go,
        :catch,
        :hit,
        :miss,
        :aborted,
        :false_alarm,
        :correct_reject,
        :lick_latency,
        :auto_rewarded,
        :initial_image_name,
        :change_image_name]
    X = hcat([trials[:, s] for s in vars]...)
    return TimeSeries(t, vars, X)
end

flashesset = (Val{:Natural_Images_Lum_Matched_set_ophys_G_2019},
              Val{:Natural_Images_Lum_Matched_set_ophys_H_2019})
function ANB.alignlfp(session, X, stimulus::Union{flashesset...}; trail = :change,
                      trials = nothing, conditions = nothing, nonconditions = nothing,
                      zero = false)
    stimulus = string(typeof(stimulus).parameters[1])
    flash = 0.250 # s see visual behavior white paper
    iti = 0.500 # s
    is = ANB.stimulusintervals(session, stimulus)
    if trail === :change # In this case the X should only be from around change trials
        is = getchangetrials(session)
        # idxs = indexin(stimuli.trials_id, ANB.getstimuli(session, stimulus).trials_id)
        # stimuli = stimuli[idxs, :]
    end
    if !isnothing(trials)
        is = is[is.trials_id .∈ (trials.trials_id,), :]
    end
    if !isnothing(conditions)
        for c in conditions
            is = is[is[:, c] .== true, :]
        end
    end
    if !isnothing(nonconditions)
        for c in nonconditions
            is = is[is[:, c] .== false, :]
        end
    end
    if trail == :onset # The programmed change time
        onsets = is.start_time
        is = [onsets[i] .. min(onsets[i + 1], onsets[i] + flash + iti)
              for i in 1:(length(onsets) - 1)]
    elseif trail == :offset # The programmed off time
        offsets = is.stop_time
        onsets = is.start_time[2:end]
        is = [offsets[i] .. min(offsets[i + 1], offsets[i] + iti)
              for i in 1:(length(offsets) - 1)]
    elseif trail == :change # The true change time, including display delay
        onsets = is.change_time_with_display_delay
        onsets = [onsets; Inf]
        is = [onsets[i] .. min(onsets[i + 1], onsets[i] + flash + iti)
              for i in 1:(length(onsets) - 1)]
    else
        is = is.interval
    end
    X = ANB.rectifytime(X)
    _X = [X[𝑡(g)] for g in is]
    _X = _X[.!isempty.(_X)]
    _X = [x for x in _X if size(x, Ti) > 100] # Remove short trials
    _X = [x[1:minimum(size.(_X, Ti)), :] for x in _X] # Catch any that are one sample too long
    _X = ANB.rectifytime.(_X; zero)
    @assert all(all(size.(x, 1) .== size(x[1], 1)) for x in _X)
    return _X
end

end
export VisualBehavior
getprobefiles(S::AbstractNWBSession; dataset = VisualBehavior) = dataset.getprobefiles(S)

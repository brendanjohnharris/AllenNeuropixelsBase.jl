using IntervalSets

export ecephyscache, getsessiontable, getprobes, getchannels, listprobes, getsessiondata, AbstractSession, Session, getid, getprobes, getfile, getprobeids, getchannels, getprobecoordinates, getstructureacronyms, getstructureids, getprobestructures, getprobe, getunits, getepochs, getstimuli, getstimulustimes

function ecephyscache()
    ecephys_project_cache.EcephysProjectCache.from_warehouse(manifest=ecephysmanifest)
end

"""
# getsessiontable — Read and return the session table from the EcephysProjectCache
    `getsessiontable()`

Read the session table data from the `EcephysProjectCache` object returned by `ecephyscache` and returns it as a `DataFrame`.
"""
function getsessiontable()
    @info "Please wait, this can take a few seconds"
    CSV.read(IOBuffer(ecephyscache().get_session_table().to_csv()), DataFrame);
end
export getsessiontable

"""
    `getprobes()`

Read the probe data from the `EcephysProjectCache` object returned by `ecephyscache()` and return it as a `DataFrame`.

## Returns
A `DataFrame` containing the probe data.
"""
function getprobes()
    CSV.read(IOBuffer(ecephyscache().get_probes().to_csv()), DataFrame);
end
export getprobes

function getchannels()
    CSV.read(IOBuffer(ecephyscache().get_channels().to_csv()), DataFrame);
end
export getchannels

listprobes(session) = getchannels.((session,), getprobeids(session))
export listprobes


function getsessiondata(session_id::Int; filter_by_validity=true, amplitude_cutoff_maximum = 0.1, presence_ratio_minimum = 0.9, isi_violations_maximum = 0.5)
    ecephyscache().get_session_data(session_id; filter_by_validity=filter_by_validity,
                            amplitude_cutoff_maximum=amplitude_cutoff_maximum,
                            presence_ratio_minimum=presence_ratio_minimum,
                            isi_violations_maximum=isi_violations_maximum)
end
export getsessiondata

abstract type AbstractSession end
export AbstractSession

struct Session <: AbstractSession
    pyObject
end
export Session
function Session(session_id::Int; kwargs...)
    Session(getsessiondata(session_id; kwargs...))
end
Session(; params...) = Session(params[:sessionid]);
getid(S::AbstractSession) = pyconvert(Int, S.pyObject.ecephys_session_id)
getprobes(S::AbstractSession) = py2df(S.pyObject.probes)
getfile(S::AbstractSession) = datadir*"Ecephys/session_"*string(getid(S))*"/session_"*string(getid(S))*".nwb"
getprobeids(S::AbstractSession) = getprobes(S)[!, :id]
function getchannels(S::AbstractSession)
    df = py2df(S.pyObject.channels)
    df.structure_acronym = df.ecephys_structure_acronym
    return df
end
function getchannels(S::AbstractSession, probeid)
    c = getchannels(S)
    c = subset(c, :probe_id=>ByRow(==(probeid)))
end
function getprobecoordinates(S::AbstractSession)
    c = subset(getchannels(S),              :anterior_posterior_ccf_coordinate => ByRow(!ismissing),
                                            :dorsal_ventral_ccf_coordinate => ByRow(!ismissing),
                                            :left_right_ccf_coordinate => ByRow(!ismissing))
    x = c[!, :anterior_posterior_ccf_coordinate]
    y = c[!, :dorsal_ventral_ccf_coordinate]
    z = c[!, :left_right_ccf_coordinate]
    return (x, y, z)
end
function getprobecoordinates(S::AbstractSession, probeid)
    c = subset(getchannels(S, probeid),              :anterior_posterior_ccf_coordinate => ByRow(!ismissing),
                                            :dorsal_ventral_ccf_coordinate => ByRow(!ismissing),
                                            :left_right_ccf_coordinate => ByRow(!ismissing))
    x = c[!, :anterior_posterior_ccf_coordinate]
    y = c[!, :dorsal_ventral_ccf_coordinate]
    z = c[!, :left_right_ccf_coordinate]
    return (x, y, z)
end

notemptyfirst(x) = length(x) > 0 ? x[1] : missing

function getstructureacronyms(channelids::Vector{Int})
    channels = getchannels()
    acronyms = Vector{Any}(undef, size(channelids))
    [acronyms[i] = notemptyfirst(channels[channels.id.==channelids[i], :structure_acronym]) for i ∈ 1:length(channelids)]
    return acronyms
end

function getstructureids(channelids::Vector{Int})
    channels = getchannels()
    acronyms = Vector{Any}(undef, size(channelids))
    [acronyms[i] = notemptyfirst(channels[channels.id.==channelids[i], :ecephys_structure_id]) for i ∈ 1:length(channelids)]
    return acronyms
end

function getprobestructures(S::AbstractSession)
    df = listprobes(S)
    acronyms = [unique(string.(d.structure_acronym)) for d in df]
    # acronyms = [a[.!ismissing.(a)] for a in acronyms]
    probeids = [unique(d.probe_id) |> first for d in df]
    return Dict(probeids .=> acronyms)
end

function getprobestructures(S::AbstractSession, structures::AbstractVector)
    D = getprobestructures(S)
    filter!(d->any(structures .∈ (last(d),)), D)
    D = Dict(k=>first(v[v.∈[structures]]) for (k, v) in D)
end

function getprobe(S::AbstractSession, structure::AbstractString)
    D = getprobestructures(S, [structure])
    return first(keys(D))
end


function getunits(; filter_by_validity=true, amplitude_cutoff_maximum = 0.1, presence_ratio_minimum = 0.9, isi_violations_maximum = 0.5)
    str = ecephyscache().get_units(filter_by_validity=filter_by_validity,
                            amplitude_cutoff_maximum=amplitude_cutoff_maximum,
                            presence_ratio_minimum=presence_ratio_minimum,
                            isi_violations_maximum=isi_violations_maximum).to_csv()
    CSV.read(IOBuffer(str), DataFrame);
end
export getunits

function getunits(structure::String; kwargs...)
    units = getunits(; kwargs...)
    units = subset(units, :structure_acronym, structure)
end
function getsessionunits(session::AbstractSession) # Much faster than the option below
    units = pyObject(session.pyObject.units.to_dataframe())
    units.structure_acronym = units.ecephys_structure_acronym
end
function getsessionunits(session::AbstractSession, structure::String)
    subset(getsessionunits(session), :structure_acronym, structure)
end
function getunits(session::AbstractSession; kwargs...)
    units = getunits(; kwargs...)
    units = subset(units, :ecephys_session_id, getid(session))
end
getunits(session::AbstractSession, structure::String) = subset(getunits(session), :structure_acronym, structure)



function getunitanalysismetricsbysessiontype(session_type; filter_by_validity=true, amplitude_cutoff_maximum = 0.1, presence_ratio_minimum = 0.9, isi_violations_maximum = 0.5) # Yeah thats python
    str = ecephyscache().get_unit_analysis_metrics_by_session_type(session_type,
                            filter_by_validity=filter_by_validity,
                            amplitude_cutoff_maximum=amplitude_cutoff_maximum,
                            presence_ratio_minimum=presence_ratio_minimum,
                            isi_violations_maximum=isi_violations_maximum).to_csv()
    CSV.read(IOBuffer(str), DataFrame);
end
export getunitanalysismetricsbysessiontype



function getallunitmetrics() # This one is really slow
    metrics1 = get_unit_analysis_metrics_by_session_type("brain_observatory_1.1",
                            amplitude_cutoff_maximum = Inf,
                            presence_ratio_minimum = -Inf,
                            isi_violations_maximum = Inf)

    metrics2 = get_unit_analysis_metrics_by_session_type("functional_connectivity",
                            amplitude_cutoff_maximum = Inf,
                            presence_ratio_minimum = -Inf,
                            isi_violations_maximum = Inf)

    vcat(analysis_metrics1, analysis_metrics2)
end
export getallunitmetrics

function getunitanalysismetrics(session::AbstractSession; annotate=true, filter_by_validity=true, kwargs...)
    str = ecephyscache().get_unit_analysis_metrics_for_session(getid(session); annotate, filter_by_validity, kwargs...)
    PyPandasDataFrame(str) |> DataFrame
end


function getstimuli(S::Session)
    str =  S.pyObject.stimulus_presentations
    PyPandasDataFrame(str) |> DataFrame
end

function getunitmetrics(session::AbstractSession)
    str = session.pyObject.units
    PyPandasDataFrame(str) |> DataFrame
end

function getstimulusname(session::AbstractSession, time::Number; stimulus_table=getstimuli(session))
    idx = findlast(stimulus_table.start_time .< time)
    if isnothing(idx)
        "blank"
    else
        stimulus_table.stimulus_name[idx]
    end
end
getstimulusname(session::AbstractSession, times; stimulus_table=getstimuli(session), kwargs...) = getstimulusname.([session], times; stimulus_table, kwargs...)


function getstimuli(S::Session, stimulusname::String)
    stimulus_table = getstimuli(S)
    df = subset(stimulus_table, :stimulus_name=>ByRow(==(stimulusname)))
end

function getstimuli(session::Session, times::Union{Tuple, UnitRange, LinRange, Vector})
    stimuli = getstimuli(session)
    idxs = [findfirst(time .< stimuli.stop_time) for time ∈ times] # Find first frame that ends after each time point
    return stimuli[idxs, :]
end

function getstimuli(session::Session, times::Interval)
    stimuli = getstimuli(session)
    start = stimuli.start_time
    stop = stimuli.stop_time
    i1 = min(findfirst(start .> minimum(times)), findfirst(stop .> minimum(times)))
    i2 = max(findlast(start .< maximum(times)), findlast(stop .< maximum(times)))
    stimuli = stimuli[i1:i2, :]
end

function getepochs(S::Session)
    p = S.pyObject.get_stimulus_epochs() # Why is this so slow
    PyPandasDataFrame(p) |> DataFrame
end

function getepochs(S::AbstractSession, stimulusname)
    epoch_table = getepochs(S)
    df = subset(epoch_table, :stimulus_name=>ByRow(==(stimulusname)))
end

function getstimulustimes(S::Session, stimulusname)
    E = getepochs(S, stimulusname)
    ts = [E.start_time E.stop_time]
    times = [x..y for (x, y) in eachrow(ts)]
end

function getstimulustimes(; params...)
    S = Session(params[:sessionid])
    getstimulustimes(S, params[:stimulus])[params[:epoch]]
end

# * Stimulus analysis

ReceptiveFieldMapping(S::AbstractSession) = stimulusmapping.ReceptiveFieldMapping(S.pyObject)

function getreceptivefield(S::AbstractSession, channel::Number)
    rfm = ReceptiveFieldMapping(S)
    rfm.get_receptive_field(channel)
end

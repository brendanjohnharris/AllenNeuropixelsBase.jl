using IntervalSets
using HDF5
using Statistics

export LFPVector, LFPMatrix, PSDMatrix, PSDVector, LogPSDVector, duration, samplingperiod, getlfp, getlfptimes, getlfpchannels, samplingrate, WaveletMatrix, LogWaveletMatrix, formatlfp, getchannels, getchanneldepths, waveletmatrix, getunitdepths, getdim, gettimes, sortbydepth, rectifytime, stimulusepochs, stimulusintervals, gaborintervals, alignlfp, logwaveletmatrix, matchlfp, joinlfp, catlfp, channels2depths, selectepochs, getchannellayers, getstreamlines

LFPVector = AbstractDimArray{T,1,Tuple{A},B} where {T,A<:DimensionalData.TimeDim,B}
LFPMatrix = AbstractDimArray{T,2,Tuple{A,B}} where {T,A<:DimensionalData.TimeDim,B<:Dim{:channel}}
export LFPMatrix, LFPVector # For simpler dispatch
dimmatrix(a::Symbol, b::Symbol) = AbstractDimArray{T,2,Tuple{A,B}} where {T,A<:Dim{a},B<:Dim{b}}
dimmatrix(a, b::Symbol) = AbstractDimArray{T,2,Tuple{A,B}} where {T,A<:a,B<:Dim{b}}
dimmatrix(a, b) = AbstractDimArray{T,2,Tuple{A,B}} where {T,A<:a,B<:b}
export dimmatrix

PSDMatrix = dimmatrix(:frequency, :channel)
PSDVector = AbstractDimArray{T,1,Tuple{A},B} where {T,A<:Dim{:frequency},B}
LogPSDVector = AbstractDimArray{T,1,Tuple{A},B} where {T,A<:Dim{:logfrequency},B}

duration(X::AbstractDimArray) = diff(extrema(dims(X, Ti)) |> collect) |> first
function samplingperiod(X::AbstractDimArray)
    if dims(X, Ti).val.data isa AbstractRange
        step(dims(X, Ti))
    else
        mean(diff(collect(dims(X, Ti))))
    end
end
samplingrate(X::AbstractDimArray) = 1 / samplingperiod(X)


WaveletMatrix = dimmatrix(Ti, :frequency) # Type for DimArrays containing wavelet transform info
LogWaveletMatrix = dimmatrix(Ti, :logfrequency) # Type for DimArrays containing wavelet transform info
export WaveletMatrix, LogWaveletMatrix
function Base.convert(::Type{LogWaveletMatrix}, x::WaveletMatrix)
    x = DimArray(x, (dims(x, Ti), Dim{:logfrequency}(log10.(dims(x, :frequency)))); metadata=metadata(x), refdims=refdims(x))
    x = x[:, .!isinf.(dims(x, :logfrequency))]
end
function Base.convert(::Type{LogPSDVector}, x::PSDVector)
    x = DimArray(x, (Dim{:logfrequency}(log10.(dims(x, :frequency))),); metadata=metadata(x), refdims=refdims(x))
    x = x[.!isinf.(dims(x, :logfrequency))]
end
function Base.convert(::Type{WaveletMatrix}, x::LogWaveletMatrix)
    x = DimArray(x, (dims(x, Ti), Dim{:frequency}(exp10.(dims(x, :logfrequency)))); metadata=metadata(x), refdims=refdims(x))
end
waveletmatrix(res::LogWaveletMatrix) = convert(WaveletMatrix, res)
logwaveletmatrix(res::WaveletMatrix) = convert(LogWaveletMatrix, res)
logpsdvector(res::PSDVector) = convert(LogPSDVector, res)











function downloadlfp(S::AbstractSession, probeid::Int)
    @assert(any(getprobeids(S) .== probeid), "Probe $probeid does not belong to session $(getid(S))")
    @assert(subset(getprobes(S), :id => ByRow(==(probeid)))[!, :has_lfp_data][1], @error "Probe $probeid does not have LFP data")
    _ = S.pyObject.get_lfp(probeid)
    return nothing
end


function structure2probe(S::AbstractSession, structure::String)
    channels = getchannels(S)
    channels = channels[.!ismissing.(channels.structure_acronym), :]
    channels = subset(channels, :structure_acronym, structure)
    @assert length(unique(channels.probe_id)) == 1
    return channels.probe_id |> unique |> first
end

function getlfptimes(session::AbstractSession, probeid)
    path = getlfppath(session, probeid)
    if !isfile(path)
        @error "Probe lfp data file does not exist. This may be becuase you have not downloaded the probe data. Try using `getlfp`"
    end

    f = h5open(path)
    times = f["acquisition"][splitext(basename(path))[1]][splitext(basename(path))[1]*"_data"]["timestamps"][:]
    close(f)
    return times
end

function getlfptimes(session::AbstractSession, probeid, idxs)
    getlfptimes(session, probeid)[idxs]
end
function getlfptimes(session::AbstractSession, probeid::Int, i::Interval)
    ts = getlfptimes(session::AbstractSession, probeid::Int)
    ts = ts[ts.∈(i,)]
end

function getlfpchannels(session::AbstractSession, probeid)
    if !isfile(getlfppath(session, probeid))
        downloadlfp(session, probeid)
    end
    f = h5open(getlfppath(session, probeid))
    channels = f["general"]["extracellular_ephys"]["electrodes"]["id"][:]
    close(f)
    return channels
end

function getprobeobj(S::AbstractSession, probeid)
    probes = S.pyObject.probes |> py2df
    probeids = probes.id
    thisprobe = findfirst(probeids .== probeid)
    probe = probes[thisprobe, :]
end

function resolvenwblfp(S::AbstractSession, probeid)
    probe = getprobeobj(S, probeid)
    res = "probe_$(probe.id)_lfp"
end

"""
Get the lfp data for a probe, providing *indices* for channels and times. See function below for indexing by channel ids and time values/intervals
"""
function _getlfp(session::AbstractSession, probeid::Int; channelidxs=1:length(getlfpchannels(session, probeid)), timeidxs=1:length(getlfptimes(session, probeid)))
    @assert(any(getprobeids(session) .== probeid), "Probe $probeid does not belong to session $(getid(session))")
    @assert(subset(getprobes(session), :id => ByRow(==(probeid)))[!, :has_lfp_data][1], @error "Probe $probeid does not have LFP data")
    path = getlfppath(session, probeid)
    if !isfile(path)
        downloadlfp(session, probeid)
    end
    if !isfile(path)
        @error "LFP data did not download correctly"
    end

    f = h5open(path)

    timedata = getlfptimes(session, probeid, timeidxs)
    dopermute = true
    channelids = getlfpchannels(session, probeid)
    channelids = channelids[channelidxs]
    res = resolvenwblfp(session, probeid)
    if (channelidxs isa Union{Int64,AbstractRange{Int64}}) & (timeidxs isa Union{Int64,AbstractRange{Int64}}) # Can use HDF5 slicing
        lfp = f["acquisition"][res][res*"_data"]["data"][channelidxs, timeidxs]
    elseif timeidxs isa Union{Int64,AbstractRange{Int64}}
        lfp = [f["acquisition"][res][res*"_data"]["data"][i, timeidxs] for i ∈ channelidxs]
        lfp = hcat(lfp...)
        dopermute = false
    else
        lfp = read(f["acquisition"][res][res*"_data"]["data"])
        lfp = lfp[channelidxs, timeidxs]
    end
    if lfp isa Vector
        lfp = reshape(lfp, 1, length(lfp))
    end
    if channelids isa Number
        channelids = [channelids]
    end
    if dopermute
        lfp = permutedims(lfp, reverse(1:ndims(lfp)))
    end
    X = DimArray(lfp, (Ti(timedata), Dim{:channel}(channelids)); metadata=Dict(:sessionid => getid(session), :probeid => probeid))
    close(f)
    return X
end


function isinvalidtime(session::AbstractSession, probeids=getprobeids(session), times=NaN)
    if hasproperty(session.pyObject, :get_invalid_times) && isempty(session.pyObject.get_invalid_times()) # No invalid times in this session!
        return false
    end
    intervals = [session.pyObject.get_invalid_times().start_time.values, session.pyObject.get_invalid_times().stop_time.values]
    intervals = [[pyconvert(Float64, intervals[1][i]), pyconvert(Float64, intervals[2][i])] for i ∈ 0:(length(intervals[1])-1)]
    tags = session.pyObject.get_invalid_times().tags.values
    # tags = vcat([pyconvert(Array{Int64}, i) for i ∈ tags]...)
    # badprobes = tags[.!isnothing.(tags)]
    badprobes = []
    if times isa Interval
        isininterval = [any(i .∈ (times,)) for i ∈ intervals]
    else
        isininterval = [any((times .> i[1]) .& (times .< i[2])) for i ∈ intervals]
    end
    return any(probeids .∈ (badprobes,)) & any(isininterval)
end
function isinvalidtime(session::AbstractNWBSession, probeids=getprobeids(session), times=NaN)
    return !pyconvert(Bool, getfile(session).invalid_times == @py None)
end

"""
This is the one you should be using. Get lfp data by channel id and time intervals or vector. Also, throw error if you try to access an invalid time interval.
"""
function getlfp(session::AbstractSession, probeid::Int; channels=getlfpchannels(session, probeid), times=ClosedInterval(extrema(getlfptimes(session, probeid))...), inbrain=false)
    if isinvalidtime(session, probeid, times)
        @error "Requested LFP data contains an invalid time..."
    end
    if inbrain isa Symbol || inbrain isa Real || inbrain
        depths = getchanneldepths(session, probeid, channels)
        if inbrain isa Real # A depth cutoff
            channels = channels[depths.>inbrain]
        elseif inbrain isa Symbol # A mode
        else # Just cutoff at the surface
            channels = channels[depths.>0]
        end
    end

    channelidxs = getlfpchannels(session, probeid)
    channelidxs = indexin(channels, channelidxs)
    channelidxs = filter(!isnothing, channelidxs)
    @assert length(channels) == length(channelidxs)

    timeidxs = getlfptimes(session, probeid)
    if !(times isa Interval) && length(times) == 2
        times = ClosedInterval(times...)
    end
    if times isa Interval
        timeidxs = findall(timeidxs .∈ (times,))
        if isempty(timeidxs)
            @warn "Session $(getid(session)), probe $probeid has no LFP for $(times)"
        end
    else
        timeidxs = indexin(times, timeidxs)
        timeidxs = filter(!isnothing, timeidxs)
        @assert length(timeidxs) == length(times)
    end

    # See if we can convert to unitranges for faster HDF5 reading via slicing
    if collect(UnitRange(extrema(timeidxs)...)) == timeidxs
        timeidxs = UnitRange(extrema(timeidxs)...)
    end
    if collect(UnitRange(extrema(channelidxs)...)) == channelidxs
        channelidxs = UnitRange(extrema(channelidxs)...)
    end

    @info "Accessing LFP data for session $(getid(session))"
    _getlfp(session, probeid; channelidxs, timeidxs)
end

"""
If you want to downsample the LFP data, its quicker to use this function and then perform slicing afterwards (since getlfp() has to check all of the time coordinates you supply, which can be slow).
"""
function getdownsampledlfp(session, probeid; downsample=100, timerange=ClosedInterval(extrema(getlfptimes(session, probeid))...), channels=getlfpchannels(session, probeid))
    if !(timerange isa Interval) && length(timerange) == 2
        timerange = ClosedInterval(timerange...)
    end
    timevals = getlfptimes(session, probeid)
    tidxs = timevals .∈ (timerange,)
    times = findfirst(tidxs):downsample:findlast(tidxs)
    _getlfp(session, probeid; timeidxs=times)[:, At(channels)]
end


"""
Now we can overload `getlfp()` to index by structure
"""
function getlfp(session::AbstractSession, probeid::Int, structures::Union{Vector{<:AbstractString},AbstractString}; kwargs...)
    if structures isa String
        structures = [structures]
    end
    channels = subset(getchannels(session, probeid), :structure_acronym => ByRow(in(structures)), skipmissing=true)
    channels = channels.id ∩ getlfpchannels(session, probeid)
    isempty(channels) && @error "No matching channels found for structure(s) $structures. Perhaps you have entered the wrong probe id?"
    getlfp(session, probeid; channels, kwargs...)
end

function getlfp(S::AbstractSession, structure::AbstractString; kwargs...)
    probeid = structure2probe(S, structure)
    getlfp(S, probeid, structure; kwargs...)
end

function getlfp(session, probeids::Vector{Int}, args...; kwargs...)
    LFP = [getlfp(session, probeid, args...; kwargs...) for probeid ∈ probeids]
end

function selectepochs(session, stimulus, epoch)
    epochs = getepochs(session, stimulus)
    isnothing(epoch) && (return epochs)
    isempty(epochs) && (@error "No '$stimulus' epoch found in session $(getid(session))")
    if !(epoch isa Symbol) && length(epoch) > 1
        qualifier = epoch[2]
        epoch = epoch[1]
        _epochs = subset(epochs, qualifier => ByRow(==(true)))
        if !isempty(_epochs)
            epochs = _epochs
        else
            error("Invalid qualifier '$qualifier' for epoch '$epoch' in this session")
        end
    end
    if epoch == :longest
        _, epoch = findmax(epochs.duration)
    elseif epoch == :first
        epoch = 1
    elseif epoch == :last
        epoch = size(epochs, 1)
    end
    epoch = epochs[epoch, :] |> DataFrame
end

function formatlfp(session::AbstractSession; probeid=nothing, tol=6, stimulus="gabors", structure="VISp", epoch=:longest, kwargs...)
    if isnothing(probeid)
        probeid = getprobe(session, structure)
    end
    if isnothing(structure)
        structure = getstructureacronyms(session, getchannels(session, probeid).id) |> unique |> skipmissing |> collect |> Vector{String}
        structure = structure[structure.!=["root"]]
    end
    if stimulus == "all"
        X = rectifytime(getlfp(session, probeid, structure; inbrain=0); tol)
    else
        epoch = selectepochs(session, stimulus, epoch)
        times = epoch.start_time .. epoch.stop_time
        X = rectifytime(getlfp(session, probeid, structure; inbrain=0, times); tol)
    end
    X = sortbydepth(session, probeid, X)
    X = DimArray(X; metadata=Dict(:sessionid => getid(session), :probeid => probeid, :stimulus => stimulus, :structure => structure))
end

function formatlfp(; sessionid=757216464, probeid=nothing, kwargs...)
    if sessionid < 1000000000
        session = Session(sessionid)
    else
        session = VisualBehavior.Session(sessionid)
    end
    formatlfp(session; probeid, kwargs...)
end
export formatlfp




function getstreamlines()
    streamlinepath = abspath(referencespacemanifest, "../laplacian_10.nrrd")
    if !isfile(streamlinepath)
        # Download it
        @info "Downloading streamline data, this may take a few minutes"
        download("https://www.dropbox.com/sh/7me5sdmyt5wcxwu/AACFY9PQ6c79AiTsP8naYZUoa/laplacian_10.nrrd?dl=1", streamlinepath)

    end
    streamlines = load(streamlinepath)
    streamlines = DimArray(streamlines.data,
        (
            Dim{:L}(getproperty(streamlines.axes[1], :val)),
            Dim{:P}(getproperty(streamlines.axes[2], :val)),
            Dim{:S}(getproperty(streamlines.axes[3], :val)),
        ))
    return streamlines
end



function getchannellayers(session, channels, cdf=getchannels(session))
    function get_layer_name(acronym)
        try
            if !contains(acronym, "VIS")
                return 0
            end
            layer = parse(Int, match(r"\d+", acronym).match)
            if layer == 3
                layer = 0
            end
            return layer
        catch
            return 0
        end
    end

    function get_structure_ids(df, annotations)
        x = floor.(Int, df.anterior_posterior_ccf_coordinate / 10)
        y = floor.(Int, df.dorsal_ventral_ccf_coordinate / 10)
        z = floor.(Int, df.left_right_ccf_coordinate / 10)

        x[x.<0] .= 0
        y[y.<0] .= 0
        z[z.<0] .= 0

        structure_ids = [annotations[_x.+1, _y.+1, _z.+1] for (_x, _y, _z) in zip(x, y, z)] # annotation volume is 1-indexed

        return structure_ids .|> Int |> vec
    end

    df = cdf[indexin(channels, cdf.id), :]
    annotations, _ = getannotationvolume(; resolution=10)
    df = df[df.anterior_posterior_ccf_coordinate.>0, :]
    x = floor.(Int, df.anterior_posterior_ccf_coordinate / 10)
    y = floor.(Int, df.dorsal_ventral_ccf_coordinate / 10)
    z = floor.(Int, df.left_right_ccf_coordinate / 10)

    structure_ids = get_structure_ids(df, annotations)
    structure_tree = Dict(v => k for (k, v) in getstructureidmap())
    structure_acronyms = [s == 0 ? "root" : getindex(structure_tree, s) for s in structure_ids]
    layers = [get_layer_name(acronym) for acronym in structure_acronyms]
    return layers, structure_acronyms
end


function getchannels(data::AbstractDimArray)
    dims(data, :channel).val
end

"""
At the moment this is just a proxy: distance along the probe to the cortical surface
"""
function getchanneldepths(session, probeid, channels; kwargs...)
    cdf = getchannels(session, probeid)
    return _getchanneldepths(cdf, channels; kwargs...)
end
function getchanneldepths(session, channels::Union{AbstractVector,Tuple}; kwargs...)
    cdf = getchannels(session) # Slightly slower than the above
    cdf = cdf[indexin(channels, cdf.id), :]
    cdfs = groupby(cdf, :probe_id)
    depths = vcat([_getchanneldepths(c, c.id; kwargs...) for c ∈ cdfs]...)
    depths = depths[indexin(channels, vcat(cdfs...).id)]
end
function _getchanneldepths(cdf, channels; method=:streamlines)
    # surfaceposition = minimum(subset(cdf, :structure_acronym=>ByRow(ismissing)).probe_vertical_position) # Minimum because tip is at 0

    if method === :dorsal_ventral
        if any(ismissing.(cdf.structure_acronym))
            surfaceposition = maximum(subset(cdf, :structure_acronym => ByRow(ismissing)).dorsal_ventral_ccf_coordinate)
        elseif any(skipmissing(cdf.structure_acronym) .== ["root"]) # For VBN files, "root" rather than "missing"
            surfaceposition = maximum(subset(cdf, :structure_acronym => ByRow(==("root"))).dorsal_ventral_ccf_coordinate)
        end
        # Assume the first `missing` channel corresponds to the surfaceprobe_vertical_position
        idxs = indexin(channels, cdf.id)[:]
        # alldepths = surfaceposition .- cdf.probe_vertical_position # in μm
        alldepths = cdf.dorsal_ventral_ccf_coordinate .- surfaceposition # in μm
        depths = fill(NaN, size(idxs))
        depths[.!isnothing.(idxs)] = alldepths[idxs[.!isnothing.(idxs)]]
    elseif method === :streamlines # This one only really works for the cortex. Anythign outside the cortex is a guess based on linear extrapolation.
        # Also note that this gives what seems to be a proportion in the cortex
        # https://www.dropbox.com/sh/7me5sdmyt5wcxwu/AACB2idSLV-F_QOG894NnZS2a?dl=0&preview=layer_mapping_example.py # The code used in the siegle 2021 paper

        streamlines = getstreamlines()

        df = cdf[cdf.anterior_posterior_ccf_coordinate.>0, :]
        df = df[indexin(channels, df.id), :]
        x = df.anterior_posterior_ccf_coordinate
        y = df.dorsal_ventral_ccf_coordinate
        z = df.left_right_ccf_coordinate


        cortical_depth = [streamlines[Dim{:L}(Near(_x)), Dim{:P}(Near(_y)), Dim{:S}(Near(_z))] for (_x, _y, _z) ∈ zip(x, y, z)] # 1-based indexing

        # Linearly extrapolate the zero depths
        _xs = findall(cortical_depth .> 0)
        _ys = cortical_depth[_xs]
        depthf(x) = first([1 x] * ([ones(length(_xs)) _xs] \ _ys))
        cortical_depth[cortical_depth.==0] .= depthf.(findall(cortical_depth .== 0))
        # f = Figure()
        # ax = Axis3(f[1, 1]; aspect=:data)
        # volume!(ax, dims(streamlines, Dim{:L})[1:100:end] |> collect,
        #     dims(streamlines, Dim{:P})[1:10:end] |> collect,
        #     dims(streamlines, Dim{:S})[1:100:end] |> collect,
        #     streamlines.data[1:100:end, 1:100:end, 1:100:end])
        # meshscatter!(ax, x, y, z, markersize=100, color=AN.getchanneldepths(session, channels; method=:dorsal_ventral))
        # current_figure()


        df.cortical_depth .= 0.0
        df[df.anterior_posterior_ccf_coordinate.>0, :cortical_depth] .= cortical_depth
        depths = df.cortical_depth

    end
    return depths
end
function getchanneldepths(session, probeid, X::LFPMatrix; kwargs...)
    channels = dims(X, Dim{:channel}) |> collect
    getchanneldepths(session, probeid, channels; kwargs...)
end
function getchanneldepths(X::LFPMatrix; kwargs...)
    @assert all(haskey.((metadata(X),), (:sessionid, :probeid)))
    S = Session(metadata(X)[:sessionid])
    getchanneldepths(S, metadata(X)[:probeid], X; kwargs...)
end

function channels2depths(session, probeid::Integer, X::AbstractDimArray, d::Integer; kwargs...)
    Y = deepcopy(X)
    _d = d
    c = dims(Y, _d) |> collect
    depths = getchanneldepths(session, probeid, c; kwargs...)
    Y = set(Y, dims(Y, _d) => Dim{:depth}(depths))
    Y = reorder(Y, dims(Y, _d) => DimensionalData.ForwardOrdered)
    return Y
end
function channels2depths(session, probeids::Vector, X::AbstractDimArray, d; kwargs...)
    Y = deepcopy(X)
    d = [d...]
    for (i, _d) in d
        probeid = probeids[i]
        c = dims(Y, _d) |> collect
        depths = getchanneldepths(session, probeid, c; kwargs...)
        Y = set(Y, dims(Y, _d) => Dim{:depth}(depths))
    end
    return Y
end

function getunitdepths(session, probeid, units; kwargs...)
    metrics = getunitanalysismetrics(session)
    check = all(units .∈ (metrics.ecephys_unit_id,))
    if !check
        @error "Some units do not have corresponding metrics"
    end
    metrics = subset(metrics, :ecephys_unit_id, units)
    channels = metrics.ecephys_channel_id
    getchanneldepths(session, probeid, channels; kwargs...)
end

getdim(X::AbstractDimArray, dim) = dims(X, dim).val
gettimes(X::AbstractDimArray) = getdim(X, Ti)


function sortbydepth(session, probeid, LFP::AbstractDimArray)
    depths = getchanneldepths(session, probeid, getchannels(LFP))
    indices = Array{Any}([1:size(LFP, i) for i in 1:length(size(LFP))])
    indices[findfirst(isa.(dims(LFP), Dim{:channel}))] = sortperm(depths)
    return LFP[indices...]
end


function rectifytime(X::AbstractDimArray; tol=6) # tol gives significant figures for rounding
    ts = gettimes(X)
    stp = ts |> diff |> mean
    err = ts |> diff |> std
    if err > stp / 10.0^(-tol)
        @warn "Time step is not approximately constant, skipping rectification"
    else
        stp = round(stp; digits=tol)
        t0, t1 = round.(extrema(ts); digits=tol)
        ts = t0:stp:t1+(10000*stp)
        ts = ts[1:size(X, Ti)] # Should be ok?
    end
    @assert length(ts) == size(X, Ti)
    X = set(X, Ti => ts)
end


function stimulusepochs(session, stim)
    stimtable = getepochs(session, stim)
    stimtable.interval = [a .. b for (a, b) in zip(stimtable.start_time, stimtable.stop_time)]
    return stimtable
end
function stimulusintervals(session, stim)
    stimtable = getstimuli(session, stim)
    stimtable.interval = [a .. b for (a, b) in zip(stimtable.start_time, stimtable.stop_time)]
    return stimtable
end

function gaborintervals(session)
    stimtable = getstimuli(session, "gabors")
    stimtable.combined_pos = sqrt.(Meta.parse.(stimtable.x_position) .^ 2 .+ Meta.parse.(stimtable.y_position) .^ 2) # Radial position of the gabor stimulus
    stimtable.interval = [a .. b for (a, b) in zip(stimtable.start_time, stimtable.stop_time)]
    return stimtable
end

function radialgaborseries(session, times)
    stimtable = getstimuli(session, "gabors")
    stimtable.combined_pos = sqrt.(Meta.parse.(stimtable.x_position) .^ 2 .+ Meta.parse.(stimtable.y_position) .^ 2) # Radial position of the gabor stimulus
    gaborseries = DimArray(zeros(length(times)), (Ti(times),))
    for pos in unique(stimtable.combined_pos)
        gabortimes = [a .. b for (a, b, x) in zip(stimtable.start_time, stimtable.stop_time, stimtable.combined_pos) if x == pos]
        for ts in gabortimes
            gaborseries[Ti(ts)] .= pos
        end
    end
    return gaborseries
end

function alignlfp(session, X, ::Val{:gabors}; x_position=nothing, y_position=nothing)
    gaborstim = gaborintervals(session)
    X = rectifytime(X)
    isnothing(x_position) || (gaborstim = gaborstim[Meta.parse.(gaborstim.x_position).==x_position, :])
    isnothing(y_position) || (gaborstim = gaborstim[Meta.parse.(gaborstim.y_position).==y_position, :])
    _X = [X[Ti(g)] for g in gaborstim.interval]
    _X = [x[1:minimum(size.(_X, Ti)), :] for x in _X] # Catch any that are one sample too long
    # _X = DimArray(mean(collect.(_X)), (Ti(step(dims(X, Ti)):step(dims(X, Ti)):step(dims(X, Ti))*minimum(size.(_X, Ti))), dims(X, Dim{:channel})))
    return _X
end

function alignlfp(session, X, ::Val{:static_gratings})
    stim = stimulusintervals(session, "static_gratings")
    stim = stim[stim.start_time.>minimum(dims(X, Ti)), :]
    stim = stim[stim.stop_time.<maximum(dims(X, Ti)), :]
    X = rectifytime(X)
    _X = [X[Ti(g)] for g in stim.interval]
    _X = [x[1:minimum(size.(_X, Ti)), :] for x in _X]
    return _X
end

"""
For flashes alignment, `trail=false` will return only the data from within the flash period. `trail=onset` will return the data from the onset of the flash to the onset of the flash through to the onset of the next flash. `trail=offset` will return the data from the offset of the flash to the onset of the next flash.
"""
function alignlfp(session, X, ::Val{:flashes}; trail=:offset)
    is = stimulusintervals(session, "flashes")
    if trail == :onset
        onsets = is.start_time
        is = [onsets[i] .. onsets[i+1] for i in 1:length(onsets)-1]
    elseif trail == :offset
        offsets = is.stop_time
        onsets = is.start_time[2:end]
        is = [offsets[i] .. onsets[i] for i in 1:length(offsets)-1]
    else
        is = is.interval
    end
    X = rectifytime(X)
    _X = [X[Ti(g)] for g in is]
    _X = [x[1:minimum(size.(_X, Ti)), :] for x in _X] # Catch any that are one sample too long
    return _X
end

# """
# For spontaneous alignment, we take each whole spontaneous interval as a trial
# """
# function alignlfp(session, X, ::Val{:spontaneous})
#     is =stimulusintervals(session, "spontaneous").interval
#     X = rectifytime(X)
#     _X = [X[Ti(g)] for g in is]
#     return _X
# end

alignlfp(session, X, stimulus::Union{String,Symbol}="gabors"; kwargs...) = alignlfp(session, X, stimulus |> Symbol |> Val; kwargs...)

"""
Adjust the times of LFP matrix Y so that they match the matrix X
"""
function matchlfp(X, Y)
    _ts = Interval(X) ∩ Interval(Y)
    ts = dims(X, Ti)
    ts = ts[last(findfirst(ts .∈ (_ts,))):last(findlast(ts .∈ (_ts,)))]
    Y = Y[Ti(Near(ts))]
    Y = DimensionalData.setdims(Y, ts)
    X = X[Ti(_ts)]
    @assert dims(X, Ti) == dims(Y, Ti)
    return (X, Y)
end

function intersectlfp(X::AbstractVector)
    Y = [rectifytime(x; tol=10) for x in X]
    ts = dims.(Y, Ti)
    ts = [Interval(extrema(t)...) for t in ts]
    int = reduce(intersect, ts)
    Y = [y[Ti(int)] for y in Y]

    ts = dims.(Y, Ti)
    s = step.(ts)
    # @assert all(s .== s[1])
    length = minimum(size.(Y, 1))
    idxs = 1:length
    Y = cat([y[idxs, :] for y in Y]..., dims=2)
end

function catlfp(X::AbstractVector)
    Y = [rectifytime(x; tol=10) for x in X]
    ts = dims.(Y, Ti)
    s = step.(ts)
    @assert all([dims(Y[1], 2)] .== dims.(Y, (2,)))
    @assert all(s .≈ s[1])
    s = s[1]
    Y = cat(Y..., dims=Ti)
    set(Y, Ti => Ti(s:s:s*size(Y, 1)))
end

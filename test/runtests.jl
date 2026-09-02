import AllenNeuropixelsBase as ANB
using AllenNeuropixelsBase.TimeseriesBase
using AllenNeuropixelsBase.NWBStream
using Test
using AllenNeuropixelsBase.DataFrames
ENV["JULIA_DEBUG"] = "CUDAExt,AllenNeuropixelsBase,AllenNeuropixels"
VCSESSIONID = 847657808
VCPROBEID = 848037574
VBSESSIONID = 1067588044

begin # * Download the files
    session_id = VCSESSIONID
    session = ANB.Session(session_id)
    ANB.getprobes(session)
    probeid = @test_nowarn ANB.getprobeids(session)[2]
    ANB.getprobestructures(session)[probeid]
    @test_nowarn ANB.listprobes(session)
    @test ANB.getepochs(session) isa DataFrame
    @test ANB.getprobes(session) isa DataFrame

    # * Check channels
    channels = ANB.getlfpchannels(session, probeid)
    cdf = ANB.getchannels(session, probeid)
    @test all(channels .∈ [cdf.id])

    LFP = ANB._getlfp(
        session, probeid;
        channelidxs = 1:length(ANB.getlfpchannels(session, probeid)),
        timeidxs = 1:length(ANB.getlfptimes(session, probeid))
    )
    LFP = []
    GC.gc() # Clean up for test runner
end

# @testset "NWBStream.jl" begin
#     f, io = @test_nowarn s3open("https://visual-behavior-neuropixels-data.s3.us-west-2.amazonaws.com/visual-behavior-neuropixels/behavior_ecephys_sessions/VCSESSIONID/ecephys_session_VCSESSIONID.nwb")
#     s3close(io)
# end

# @testset "Stream Visual Behavior" begin
#     st = @test_nowarn ANB.VisualBehavior.getsessiontable()
#     @test st isa DataFrame
#     session_id = st[2, :ecephys_session_id]
#     url = ANB.VisualBehavior.getsessionfile(session_id)
#     file, io = s3open(url)
#     ANB.behavior_ecephys_session.BehaviorEcephysSession.from_nwb(file)
#     session = ANB.VisualBehavior.S3Session(session_id)
#     ANB.initialize!(session)

#     ANB.getprobes(session)
#     ANB.getprobeids(session)
#     @test_nowarn ANB.listprobes(session)
#     @test_nowarn ANB.getepochs(session)

# end

@testset "sortbydepth keeps data with its labels" begin
    # Allen channel ids descend as cortical depth increases, so sorting by depth is a real
    # permutation of the channel axis. It must move the data and the lookup together; a previous
    # implementation moved only the lookup, mirroring every LFP against its own channel labels.
    ids = [40, 30, 20, 10]
    X = ToolsArray(Float64.(reshape(1:12, 3, 4)), (𝑡(1:3), ANB.Chan(ids)))
    ownsamples(Y) = all(parent(Y)[:, j] == parent(X)[:, findfirst(==(c), ids)]
                        for (j, c) in enumerate(collect(lookup(Y, ANB.Chan))))

    Y = ANB.sortbydepth(X, [4.0, 3.0, 2.0, 1.0])   # deepest channel first: a full reversal
    @test collect(lookup(Y, ANB.Chan)) == reverse(ids)
    @test ownsamples(Y)

    Z = ANB.sortbydepth(X, [2.0, 4.0, 1.0, 3.0])   # not a reversal, so order alone cannot pass
    @test collect(lookup(Z, ANB.Chan)) == [20, 40, 10, 30]
    @test ownsamples(Z)

    @test ANB.sortbydepth(X, [1.0, 2.0, 3.0, 4.0]) === X   # already sorted: untouched
    @test_throws ErrorException ANB.sortbydepth(X, [1.0, 2.0])
end

@testset "cached channel depths" begin
    # `addchanneldepths` writes `:depths`/`:depthmethod`; the readers used to test `:depth`/
    # `:depth_method` against an undefined `method`, so the cache was dead code.
    ids = [40, 30, 20, 10]
    X = ToolsArray(Float64.(reshape(1:12, 3, 4)), (𝑡(1:3), ANB.Chan(ids));
                   metadata = Dict(:depths => Dict(40 => 4.0, 30 => 3.0, 20 => 2.0, 10 => 1.0),
                                   :depthmethod => :probe))
    @test ANB.getchanneldepths(X; method = :probe) == [4.0, 3.0, 2.0, 1.0]

    Y = ANB.sortbydepth(X, ANB.getchanneldepths(X; method = :probe))
    @test collect(lookup(Y, ANB.Chan)) == [10, 20, 30, 40]
    @test ANB.getchanneldepths(Y; method = :probe) == [1.0, 2.0, 3.0, 4.0] # follows the lookup

    @test ANB._cacheddepths(X, :streamlines) === nothing # a different method must not reuse it
    @test ANB._cacheddepths(ToolsArray(Float64.(reshape(1:12, 3, 4)),
                                       (𝑡(1:3), ANB.Chan(ids))), :probe) === nothing
end

@testset "Format LFP" begin
    params = (;
        sessionid = VCSESSIONID,
        stimulus = "spontaneous",
        probeid = VCPROBEID,
        structure = "VISl",
        epoch = :longest,
        pass = (1, 100),
    )
    X = ANB.formatlfp(; params...)
    @test X isa ANB.LFPMatrix
    @test X isa TimeseriesBase.RegularTimeseries

    X = []
    GC.gc()
end

@testset "LFP columns match the NWB file" begin
    # The only reference that does not go through the read pipeline is the file. The
    # ElectricalSeries electrode region is the identity, so data row r belongs to the electrode
    # `getlfpchannels(...)[r]`. Asking for an ascending, contiguous run of channels is what used to
    # trigger the mirror, since it takes the range/hyperslab path and a non-trivial depth sort.
    probeid = VCPROBEID
    fileids = ANB.getlfpchannels(session, probeid)
    ts = ANB.getlfptimes(session, probeid)
    rows = 1:min(24, length(fileids))
    tidx = (length(ts) ÷ 2):(length(ts) ÷ 2 + 999)
    X = ANB.getlfp(session, probeid; channels = fileids[rows],
                   times = ts[first(tidx)] .. ts[last(tidx)])
    Z = ANB.h5open(ANB.getlfppath(session, probeid)) do f
        r = "probe_$(probeid)_lfp"
        f["acquisition"][r][r * "_data"]["data"][rows, tidx]
    end
    A = parent(X)
    m = min(size(A, 1), size(Z, 2))
    @test length(lookup(X, ANB.Chan)) == length(rows)
    @test all(A[1:m, j] ≈ Z[findfirst(==(c), fileids[rows]), 1:m]
              for (j, c) in enumerate(collect(lookup(X, ANB.Chan))))

    X = Z = A = []
    GC.gc()
end

# @testset "HybridSession" begin
#     st = @test_nowarn ANB.VisualBehavior.getsessiontable()
#     @test st isa DataFrame
#     session_id = st[2, :ecephys_session_id]
#     session = ANB.VisualBehavior.S3Session(session_id)
#     ANB.initialize!(session)

#     ANB.getprobes(session)
#     ANB.getprobeids(session)
#     @test_nowarn ANB.listprobes(session)
#     @test_nowarn ANB.getepochs(session)
# end

@testset "Visual Behavior" begin

    IS_CI = get(ENV, "CI", "false") == "true"

    if !IS_CI
        # st = @test_nowarn ANB.VisualBehavior.getsessiontable()
        # @test st isa DataFrame
        # session_id = st[end, :ecephys_session_id]
        session_id = VBSESSIONID
        session = ANB.Session(session_id)

        # test_file = "/home/brendan/OneDrive/Masters/Code/Vortices/Julia/AllenSDK/test/ecephys_session_VBSESSIONID.nwb"
        # f = ANB.behavior_ecephys_session.BehaviorSession.from_nwb_path(test_file)

        probeid = @test_nowarn ANB.getprobeids(session)[2]
        @test_nowarn ANB.getprobestructures(session)[probeid]
        @test_nowarn ANB.listprobes(session)
        @test_nowarn ANB.getepochs(session)
        @test_nowarn ANB.getprobes(session)

        # Now try to get some LFP data
        @test ANB._getlfp(
            session, probeid;
            channelidxs = 1:length(ANB.getlfpchannels(session, probeid)),
            timeidxs = 1:length(ANB.getlfptimes(session, probeid))
        ) isa
            TimeseriesBase.IrregularTimeseries

        structure = ANB.getprobestructures(session)[probeid]
        structure = structure[occursin.(("VIS",), string.(structure)) |> findfirst]

        GC.gc()

        channels = @test_nowarn ANB.getlfpchannels(session, probeid)
        cdf = @test_nowarn ANB.getchannels(session, probeid)
        ANB._getchanneldepths(cdf, channels)
        depths = @test_nowarn ANB.getchanneldepths(session, probeid, channels)

        x = ANB.getlfp(session, structure)
        @test x isa TimeseriesBase.IrregularTimeseries
        @test x isa TimeseriesBase.MultivariateTimeseries
        a = ANB.formatlfp(
            session; probeid, stimulus = "spontaneous", structure = structure,
            epoch = :longest
        )
        b = ANB.formatlfp(;
            sessionid = session_id, probeid, stimulus = "spontaneous",
            structure = structure, epoch = :longest
        ) # Slower, has to build the session
        @assert a == b

        x = a = b = []
        GC.gc()

        # Test behavior data
        S = session
        @test ANB.gettrials(S) isa DataFrame
        @test ANB.getlicks(S) isa DataFrame
        @test ANB.getrewards(S) isa DataFrame
        @test ANB.getstimuli(S) isa DataFrame
        @test ANB.geteyetracking(S) isa DataFrame
    end
end

# s3clear()

return

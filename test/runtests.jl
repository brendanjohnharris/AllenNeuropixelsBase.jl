import AllenNeuropixelsBase as ANB
using AllenNeuropixelsBase.NWBStream
using Test
using AllenNeuropixelsBase.DataFrames
ENV["JULIA_DEBUG"] = "CUDAExt,AllenNeuropixelsBase,AllenNeuropixels"

begin # * Download the files
    params = (;
        sessionid=1044385384,
        stimulus="spontaneous",
        probeid=1044506935,
        structure="VISl",
        epoch=1,
        pass=(1, 100))
    X = ANB.formatlfp(; params...)

    session_id = 1152811536
    session = ANB.Session(session_id)
    ANB.getprobes(session)
    probeid = @test_nowarn ANB.getprobeids(session)[2]
    ANB.getprobestructures(session)[probeid]
    ANB.listprobes(session)
    ANB.getepochs(session)
    ANB.getprobes(session)
    ANB._getlfp(session, probeid;
        channelidxs=1:length(ANB.getlfpchannels(session, probeid)),
        timeidxs=1:length(ANB.getlfptimes(session, probeid)))
end

# @testset "NWBStream.jl" begin
#     f, io = @test_nowarn s3open("https://visual-behavior-neuropixels-data.s3.us-west-2.amazonaws.com/visual-behavior-neuropixels/behavior_ecephys_sessions/1044385384/ecephys_session_1044385384.nwb")
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

@testset "Base" begin
    params = (;
        sessionid=1044385384,
        stimulus="spontaneous",
        probeid=1044506935,
        structure="VISl",
        epoch=1,
        pass=(1, 100))
    X = ANB.formatlfp(; params...)
    @test X isa ANB.LFPMatrix
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

import AllenNeuropixelsBase as ANB
using Test
@testset "Visual Behavior" begin
    # st = @test_nowarn ANB.VisualBehavior.getsessiontable()
    # @test st isa DataFrame
    # session_id = st[end, :ecephys_session_id]
    session_id = 1152811536
    session = ANB.Session(session_id)

    # test_file = "/home/brendan/OneDrive/Masters/Code/Vortices/Julia/AllenSDK/test/ecephys_session_1152811536.nwb"
    # f = ANB.behavior_ecephys_session.BehaviorSession.from_nwb_path(test_file)

    @test_nowarn ANB.getprobes(session)
    probeid = @test_nowarn ANB.getprobeids(session)[2]
    @test_nowarn ANB.getprobestructures(session)[probeid]
    @test_nowarn ANB.listprobes(session)
    @test_nowarn ANB.getepochs(session)
    @test_nowarn ANB.getprobes(session)

    # Now try to get some LFP data
    ANB._getlfp(session, probeid;
        channelidxs=1:length(ANB.getlfpchannels(session, probeid)),
        timeidxs=1:length(ANB.getlfptimes(session, probeid)))

    structure = ANB.getprobestructures(session)[probeid]
    structure = structure[occursin.(("VIS",), string.(structure))|>findfirst]

    channels = @test_nowarn ANB.getlfpchannels(session, probeid)
    cdf = @test_nowarn ANB.getchannels(session, probeid)
    ANB._getchanneldepths(cdf, channels)
    depths = @test_nowarn ANB.getchanneldepths(session, probeid, channels)

    ANB.getlfp(session, structure)
    a = ANB.formatlfp(session; probeid, stimulus="spontaneous", structure=structure,
        epoch=:longest)
    b = ANB.formatlfp(; sessionid=session_id, probeid, stimulus="spontaneous",
        structure=structure, epoch=:longest) # Slower, has to build the session
    @assert a == b

    # Test behavior data
    S = session
    @test ANB.gettrials(S) isa DataFrame
    @test ANB.getlicks(S) isa DataFrame
    @test ANB.getrewards(S) isa DataFrame
    @test ANB.getstimuli(S) isa DataFrame
    @test ANB.geteyetracking(S) isa DataFrame
end

# s3clear()

return

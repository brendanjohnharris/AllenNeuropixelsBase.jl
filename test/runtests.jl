import AllenNeuropixelsBase as ANB
using AllenNeuropixelsBase.NWBStream
using Test
using AllenNeuropixelsBase.DataFrames

@testset "NWBStream.jl" begin
    f, io = s3open("https://visual-behavior-neuropixels-data.s3.us-west-2.amazonaws.com/visual-behavior-neuropixels/behavior_ecephys_sessions/1044385384/ecephys_session_1044385384.nwb")
    s3close(io)
end

@testset "Stream Visual Behavior" begin
    st = @test_nowarn ANB.VisualBehavior.getsessiontable()
    @test st isa DataFrame
    session_id = st[2, :ecephys_session_id]
    url = ANB.VisualBehavior.getsessionfile(session_id)
    file, io = s3open(url)
    ANB.behavior_ecephys_session.BehaviorEcephysSession.from_nwb(file)
    session = ANB.VisualBehavior.Session(session_id)
    ANB.initialize!(session)

    ANB.getprobes(session)
    ANB.getprobeids(session)
    @test_nowarn ANB.listprobes(session)
    @test_nowarn ANB.getepochs(session)

end

# @testset "Base" begin
#     params = (;
#         sessionid = 1044385384,
#         stimulus = "gabors",
#         probeid = 769322751, # VISl # 769322749, # VISp #
#         structure = "VISp",
#         epoch = 1,
#         pass = (1, 100)
# )
#     formatlfp(; params...)
# end

@testset "HybridSession" begin
    st = @test_nowarn ANB.VisualBehavior.getsessiontable()
    @test st isa DataFrame
    session_id = st[2, :ecephys_session_id]
    session = ANB.HybridSession(session_id)
    ANB.initialize!(session)

    ANB.getprobes(session)
    ANB.getprobeids(session)
    @test_nowarn ANB.listprobes(session)
    @test_nowarn ANB.getepochs(session)
end

if false
import AllenNeuropixelsBase as ANB
using Test
@testset "Visual Behavior" begin
    st = @test_nowarn ANB.VisualBehavior.getsessiontable()
    # @test st isa DataFrame
    session_id = st[end, :ecephys_session_id]
    session = ANB.Session(session_id)

    @test_nowarn ANB.getprobes(session)
    probeid = @test_nowarn ANB.getprobeids(session)[2]
    @test_nowarn ANB.getprobestructures(session)[probeid]
    @test_nowarn ANB.listprobes(session)
    @test_nowarn ANB.getepochs(session)
    @test_nowarn ANB.getprobes(session)

    # Now try to get some LFP data
    ANB._getlfp(session, probeid; channelidxs=1:length(ANB.getlfpchannels(session, probeid)), timeidxs=1:length(getlfptimes(session, probeid)))

    ANB.getlfp(session, "VISp")
end
end

s3clear()

return

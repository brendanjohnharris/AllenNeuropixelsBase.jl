import AllenNeuropixelsBase as ANB
using AllenNeuropixelsBase.NWBS3
using Test
using AllenNeuropixelsBase.DataFrames

@testset "NWBS3.jl" begin
    f, io = s3open("https://visual-behavior-neuropixels-data.s3.us-west-2.amazonaws.com/visual-behavior-neuropixels/behavior_ecephys_sessions/1044385384/ecephys_session_1044385384.nwb")
    s3close(io)
end

@testset "Visual Behaviour" begin
    st = @test_nowarn ANB.VisualBehavior.getsessiontable()
    @test st isa DataFrame
    sessionid = st[2, :ecephys_session_id]
    session = @test_nowarn ANB.VisualBehavior.Session(sessionid)
    ANB.initialize!(session)

    @test_nowarn ANB.listprobes(session)
    @test_nowarn getepochs(session)
end

s3clear()

return

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

@testset "Base" begin
    params = (;
        sessionid = 1044385384,
        stimulus = "gabors",
        probeid = 769322751, # VISl # 769322749, # VISp #
        structure = "VISp",
        epoch = 1,
        surrogate = AP(1250),
        pass = (1, 100)
)
    formatlfp(; sessionid=757216464, probeid=769322749, stimulus="gabors", structure="VISp", epoch=1, kwargs...)
end


s3clear()

return

using NWBS3
import AllenNeuropixelsBase as ANB
using Test
using DataFrames

@testset "NWBS3.jl" begin
    s3open("https://visual-behavior-neuropixels-data.s3.us-west-2.amazonaws.com/visual-behavior-neuropixels/behavior_ecephys_sessions/1044385384/ecephys_session_1044385384.nwb")
end

# @testset "Visual Behaviour" begin
#     st = @test_nowarn ANB.VisualBehaviour.getsessiontable()
#     @test st isa DataFrame
#     sessionid = st[2, :ecephys_session_id]
#     session = @test_nowarn ANB.VisualBehaviour.Session(sessionid)
# end

import AllenNeuropixelsBase as ANB
using Test
using NWBS3
using DataFrames

@testset "Visual Behaviour" begin

st = @test_nowarn ANB.VisualBehaviour.getsessiontable()
@test st isa DataFrame
sessionid = st[2, :ecephys_session_id]
session = @test_nowarn ANB.VisualBehaviour.Session(sessionid)

end

using Colors
using GeometryBasics
using FileIO

export referencespacecache, getstructuretree, getstructurecolor, getstructurecolors, getstructuretreedepth, getstructurename, getstructuremesh, getstructureidmap

function referencespacecache(reference_space_key="annotation/ccf_2017"; resolution=25, manifest=referencespacemanifest)
    reference_space_cache.ReferenceSpaceCache(resolution, reference_space_key, manifest=manifest)
end

function getstructuretree(space=referencespacecache(), id=1) # Default id is mouse brain atlas
    id = Int(id)
    space.get_structure_tree(structure_graph_id=id)
end

function getstructurecolor(id::Union{Number,Missing}; tree=getstructuretree())
    !ismissing(id) || return RGB(0.0, 0.0, 0.0)
    id = Int(id)
    d = tree.get_structures_by_id([id])[0]
    !isnothing(d) || return RGB(0.0, 0.0, 0.0)
    trip = pyconvert(Tuple, d["rgb_triplet"])
    c = RGB((trip ./ 256)...)
end
function getstructurecolors(ids::AbstractVector)
    tree = getstructuretree()
    getstructurecolor.(ids; tree)
end

function getstructuretreedepth(id) # How many parents does this structure have (+1)?
    id = Int(id)
    tree = getstructuretree()
    d = tree.get_structures_by_id([id])[1]
    length(d["structure_id_path"])
end


function getstructurename(id::Number)
    id = Int(id)
    tree = getstructuretree()
    d = tree.get_structures_by_id([id])[1]
    d["name"]
end
function getstructurename(acronym::String)
    tree = getstructuretree()
    d = tree.get_structures_by_acronym([acronym])[1]
    d["name"]
end


function getstructureid(acronym::String)
    tree = getstructuretree()
    d = tree.get_structures_by_acronym([acronym])[1]
    d["id"]
end


function getallstructureids(args...)
    t = getstructuretree(args...)
    d = t.get_name_map()
    ids = keys(d)
end

#! Can do the rest of this dict by eval?

function buildreferencespace(tree=getstructuretree(), annotation=getannotationvolume(), resolution=(25, 25, 25))
    if annotation isa Tuple
        annotation = annotation[1]
    end
    rsp = reference_space.ReferenceSpace(tree, annotation, resolution)
end

function buildstructuremask(id, referencespace=buildreferencespace())
    referencespace.make_structure_mask([id])
end

function getstructuremeshfile(id, r=referencespacecache())
    filename = pyconvert(String, r.get_cache_path(nothing, r.STRUCTURE_MESH_KEY, r.reference_space_key, id))
    r.api.download_structure_mesh(id, r.reference_space_key, filename, strategy="lazy")
    return filename
end

function getstructureidmap()
    t = getstructuretree()
    D = t.get_id_acronym_map() |> Dict
end

function getstructuremesh(id, referencespace=referencespacecache(); hemisphere=:both)
    m = load(getstructuremeshfile(id, referencespace))

    # Only take the right hemisphere, z > midpoint
    if hemisphere === :right
        idxs = findall(last.(m.position) .> median(last.(m.position)))
    elseif hemisphere == :left
        idxs = findall(last.(m.position) .< median(last.(m.position)))
    else
        idxs = 1:length(m.position)
    end
    fs = faces(m)
    checkf(f) = all([(_f.i |> Int) ∈ idxs for _f in f]) # filter out faces that arent in the chosen hemisphere
    fs = [f for f in fs if checkf(f)]
    m = GeometryBasics.Mesh(coordinates(m), fs)
    m
end

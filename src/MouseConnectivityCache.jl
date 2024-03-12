export mouseconnectivitycache, getannotationvolume, getstructuremask, gettemplatevolume,
       getreferencespace

function mouseconnectivitycache(; kwargs...)
    mcc = mouse_connectivity_cache.MouseConnectivityCache(;
                                                          manifest_file = mouseconnectivitymanifest,
                                                          kwargs...)
end
export MouseConnectivityCache

function getannotationvolume(; resolution = 10, kwargs...)
    # vol, info = mouseconnectivitycache(; resolution, kwargs...).get_annotation_volume()
    rspc = referencespacecache(; resolution, kwargs...)
    vol, info = rspc.get_annotation_volume()
    return pyconvert(Array, vol), Dict(info)
end

function getstructuremask(structureid::Number)
    strucutreid = Int(structureid)
    mask, info = mouseconnectivitycache().get_structure_mask(structureid)
    return pyconvert(Array, mask), Dict(info)
end

function gettemplatevolume()
    template, template_info = mouseconnectivitycache().get_template_volume()
    return pyconvert(Array, template), Dict(template_info)
end

function getreferencespace()
    reference, reference_info = mouseconnectivitycache().get_reference_space()
    return pyconvert(Array, reference), Dict(reference_info)
end

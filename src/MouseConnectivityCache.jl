export mouseconnectivitycache, getannotationvolume, getstructuremask, gettemplatevolume, getreferencespace

function mouseconnectivitycache()
    mcc = mouse_connectivity_cache.MouseConnectivityCache(manifest_file=mouseconnectivitymanifest)
end
export MouseConnectivityCache


function getannotationvolume()
    vol, info = mouseconnectivitycache().get_annotation_volume()
    return pyconvert(Array, vol), Dict(info)
end

function getstructuremask(structureid::Number)
    strucutreid = Int(structureid)
    mask, info =  mouseconnectivitycache().get_structure_mask(structureid)
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

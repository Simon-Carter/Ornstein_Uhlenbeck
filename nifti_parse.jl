using NIfTI
voxel_space = niread("sample_data/processed_voxel_space/sub-108_ses-bhb_task-rest_run-1.nii");
println(size(ni))
atlas = niread("sample_data/ROIs_300inVol_MNI.nii");
println(size(atlas))

# create an array the same size as the atlas
d = findall(x -> x == 1.0, atlas.raw)

d = getindex.(d, [1 2 3 4])

#transformations

a = getaffine(atlas.header)
b = inv(getaffine(voxel_space.header))
transform = b*a

data = []

temp = Array{Float64}(undef,size(d, 1))

td = (transform * (d'))'

td = round.(Int, td[:, 1:3])

temp = Array{Float64}(undef,size(d, 1))
for i in eachindex(temp)
    temp[i] = voxel_space[td[i,:]..., 1]
end



#push!(data, Array{Float64}(undef,size(d, 1)))
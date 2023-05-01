using NIfTI
using Statistics
voxel_space = niread("sample_data/processed_voxel_space/sub-108_ses-bhb_task-rest_run-1.nii");
atlas = niread("sample_data/ROIs_300inVol_MNI.nii");
println(size(atlas))

function get_roi(voxel_space, region)
    d = findall(x -> x == region, atlas.raw)

    d = getindex.(d, [1 2 3 4])
    d[:,1:3] = d[:,1:3] .- 1

    #transformations

    a = getaffine(atlas.header)
    b = inv(getaffine(voxel_space.header))
    transform = b*a

    data = []

    temp = Array{Float64}(undef,size(d, 1))

    td = (transform * (d'))'

    td = round.(Int, td[:, 1:3])

    td = td .+ 1

    temp = Array{Float64}(undef,size(d, 1))
    k=1
    for k in 1:size(voxel_space, 4)
        for i in eachindex(temp)
            temp[i] = voxel_space[td[i,:]..., k]
        end

       #println(k, temp)

        push!(data, copy(temp))

    end

    return reduce(hcat,data)'
end

total = [get_roi(voxel_space, i) for i in 1:300];

total_dm = total[81:135]

#de-mean the data

total_deme1 = [mean(i, dims=1) for i in total_dm]
total_deme = [total_dm[i] .- total_deme1[i] for i in eachindex(total_dm)]

#seperate each region in a concatination of data
total_concat = [reshape(i, size(i,1)* size(i,2),1) for i in total_deme]

var_cocnat = var.(total_concat)

# get the average value for each ROI
total_mean = mean.(total_deme, dims=2)

var_mean = var.(total_mean)

var_diff = var_concat - var_mean
var_therm = [(size(i,2) / (size(i,2) - 1)) for i in total_deme]
var_therm = var_therm .* var_diff
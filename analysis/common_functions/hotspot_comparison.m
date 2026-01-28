% compare 2 hotspots 
% Mai-Anh Vu, 2025
% 
% NOTE: this assumes same voxel size and same voxel coordinate spacing
% hotspot1_rand = cell array of voxel indices of random hotspot1 volumes
% hotspot1_vox = indices of voxels in hotspot1
% hotspot1_DV = vector of DV coords corresponding to hotspot1 map rows
% hotspot1_ML = vector of ML coords corresponding to hotspot1 map columns
% hotspot1_AP = vector of AP coords corresponding to hotspot1 map slices

function output = hotspot_comparison(...
    hotspot1_rand,hotspot1_vox,hotspot1_DV,hotspot1_ML,hotspot1_AP,...
    hotspot2_rand,hotspot2_vox,hotspot2_DV,hotspot2_ML,hotspot2_AP)

% round 
hotspot1_DV = round(hotspot1_DV,2);
hotspot1_ML = round(hotspot1_ML,2);
hotspot1_AP = round(hotspot1_AP,2);
hotspot2_DV = round(hotspot2_DV,2);
hotspot2_ML = round(hotspot2_ML,2);
hotspot2_AP = round(hotspot2_AP,2);

% get union
DV_coords = union(hotspot1_DV,hotspot2_DV);
ML_coords = union(hotspot1_ML,hotspot2_ML);
AP_coords = union(hotspot1_AP,hotspot2_AP);
dim_union = [numel(DV_coords) numel(ML_coords) numel(AP_coords)];

% get offset
DV1_start = find(DV_coords==min(hotspot1_DV));
ML1_start = find(ML_coords==min(hotspot1_ML));
AP1_start = find(AP_coords==min(hotspot1_AP));
DV2_start = find(DV_coords==min(hotspot2_DV));
ML2_start = find(ML_coords==min(hotspot2_ML));
AP2_start = find(AP_coords==min(hotspot2_AP));

% hotspot1: recalculate voxel indices and rand voxel indices
dim1 = [numel(hotspot1_DV),numel(hotspot1_ML),numel(hotspot1_AP)];
% convert to subscripts, adjust, and then convert back to indices
[r,c,s] = ind2sub(dim1,hotspot1_vox);
r = r+DV1_start-1;
c = c+ML1_start-1;
s = s+AP1_start-1;
vox1 = sub2ind(dim_union,r,c,s);
rand1 = cell(size(hotspot1_rand));
for i = 1:numel(hotspot1_rand)
    [r,c,s] = ind2sub(dim1,hotspot1_rand{i});
    r = r+DV1_start-1;
    c = c+ML1_start-1;
    s = s+AP1_start-1;
    rand1{i} = sub2ind(dim_union,r,c,s);
end

% hotspot2: recalculate voxel indices and rand voxel indices
dim2 = [numel(hotspot2_DV),numel(hotspot2_ML),numel(hotspot2_AP)];
% convert to subscripts, adjust, and then convert back to indices
[r,c,s] = ind2sub(dim2,hotspot2_vox);
r = r+DV2_start-1;
c = c+ML2_start-1;
s = s+AP2_start-1;
vox2 = sub2ind(dim_union,r,c,s);
rand2 = cell(size(hotspot2_rand));
for i = 1:numel(hotspot2_rand)
    [r,c,s] = ind2sub(dim1,hotspot2_rand{i});
    r = r+DV2_start-1;
    c = c+ML2_start-1;
    s = s+AP2_start-1;
    rand2{i} = sub2ind(dim_union,r,c,s);
end

% now calculate overlap 
rand_overlap = arrayfun(@(x) numel(intersect(rand1{x},rand2{x})),1:numel(rand1));
actual_overlap = numel(intersect(vox1,vox2));
p_overlap = sum(rand_overlap>=actual_overlap)/numel(rand_overlap);
    
% output
output = struct;
output.n = [numel(hotspot1_vox) numel(hotspot2_vox)];
output.overlap = actual_overlap;
output.overlap_p = p_overlap;

end
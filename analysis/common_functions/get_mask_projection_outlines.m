% Mai-Anh Vu, 2025
%
% inputs
% str is output from get_striatum_vol_mask
% mask is a volume the size of str.striatum_mask
%
function output = get_mask_projection_outlines(mask,str,varargin)

%%%  parse optional inputs %%%
ip = inputParser;
% directories and atlas: see https://github.com/HoweLab/MultifiberLocalization    
ip.addParameter('apply_str_mask',1); 
ip.addParameter ('proj_orientations',{'axial','sagittal','coronal'});
ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end
if apply_str_mask==1
    mask(str.striatum_mask==0) = 0;
end

% some hard-coded info
str_coord_order.axial = [2 3];
str_coord_order.sagittal = [3 1];
str_coord_order.coronal = [1 2];

for p = 1:numel(proj_orientations)
    output.(proj_orientations{p}) = [];
    these_coords = str_coord_order.(proj_orientations{p});
    
    this_vol = get_volume_projection(mask,proj_orientations{p},...
        'mask',str.striatum_mask,'proj','max'); % max projection
    b = bwboundaries(this_vol);
    for i = 1:numel(b)
        tmp_x = str.info.(str.info.dimension_order{1,these_coords(1)});
        tmp_x = tmp_x(b{i}(:,2));
        tmp_y = str.info.(str.info.dimension_order{1,these_coords(2)});            
        tmp_y = tmp_y(b{i}(:,1));
        output.(proj_orientations{p}) = ...
            [output.(proj_orientations{p});{tmp_x} {tmp_y}];
    end
end
    
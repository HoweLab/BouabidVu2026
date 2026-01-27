% handy for masking smoothed volume maps
% Mai-Anh Vu, 2025/07/08
function output = get_striatum_vol_mask(AP_range,ML_range,DV_range,voxel_size,varargin)

    %%%  parse optional inputs %%%
    ip = inputParser;
    % directories and atlas: see https://github.com/HoweLab/MultifiberLocalization    
    ip.addParameter('atlas_dir','C:\Users\maianhvu\Documents\MATLAB\Scripts\MRIAtlas\CCF');
    ip.addParameter('factor_str_mask_into_range',1);
    ip.addParameter('voxel_size_rounding',1); % e.g., if voxel size is 0.05mm, round such that voxel boundaries are multiples of 0.05
    ip.parse(varargin{:});
    for j=fields(ip.Results)'
        eval([j{1} '=ip.Results.' j{1} ';']);
    end
    
    % get striatum mask
    [str_mask,ccf_key] = get_str_mask(atlas_dir);
    % keep right side
    str_mask(:,ccf_key.ML.ML<0,:) = 0;
    
    % include striatum boundaries if beyond input boundaries
    if factor_str_mask_into_range == 1
        strDV = sort(ccf_key.DV.DV(...
            [find(max(str_mask,[],[2 3])==1,1,'first'),...
            find(max(str_mask,[],[2 3])==1,1,'last')]));
        strAP = sort(ccf_key.AP.AP(...
            [find(max(str_mask,[],[1 2])==1,1,'first'),...
            find(max(str_mask,[],[1 2])==1,1,'last')]));
        strML = sort(ccf_key.ML.ML(...
            [find(max(str_mask,[],[1 3])==1,1,'first'),...
            find(max(str_mask,[],[1 3])==1,1,'last')]));
        % coordinate ranges based on atlas and mouse fibers
        AP_range = [min([AP_range(1) strAP(1)]) max([AP_range(2) strAP(2)])];
        ML_range = [min([ML_range(1) strML(1)]) max([ML_range(2) strML(2)])];
        DV_range = [min([DV_range(1) strDV(1)]) max([DV_range(2) strDV(2)])];
    end
    % to nearest voxel-size (e.g., to nearest 50-micron)
    if voxel_size_rounding == 1        
        DV_range(1) = floor(DV_range(1)*(1/voxel_size))/(1/voxel_size);
        DV_range(2) = ceil(DV_range(2)*(1/voxel_size))/(1/voxel_size);
        AP_range(1) = floor(AP_range(1)*(1/voxel_size))/(1/voxel_size);
        AP_range(2) = ceil(AP_range(2)*(1/voxel_size))/(1/voxel_size);
        ML_range(1) = floor(ML_range(1)*(1/voxel_size))/(1/voxel_size);
        ML_range(2) = ceil(ML_range(2)*(1/voxel_size))/(1/voxel_size);
    end

    % find the voxels in the 10micron mask corresponding to our desired
    % voxels
    vox_ML = ML_range(1):voxel_size:ML_range(2);
    vox_DV = DV_range(1):voxel_size:DV_range(2);
    vox_AP = AP_range(1):voxel_size:AP_range(2);
    [~,ML_idx] = arrayfun(@(x) min(abs(x-ccf_key.ML.ML)),vox_ML);
    [~,AP_idx] = arrayfun(@(x) min(abs(x-ccf_key.AP.AP)),vox_AP);
    [~,DV_idx] = arrayfun(@(x) min(abs(x-ccf_key.DV.DV)),vox_DV);
    [xx,yy,zz] = meshgrid(ML_idx,DV_idx,AP_idx);
    
    % output
    output = struct;
    output.info.AP = vox_AP;
    output.info.ML = vox_ML;
    output.info.DV = vox_DV;
    output.info.dimension_order = {'DV','ML','AP';'y','x','z'};
    output.info.voxel_size = voxel_size;
    output.striatum_mask = reshape(str_mask(sub2ind(size(str_mask),yy(:),xx(:),zz(:))),numel(DV_idx),numel(ML_idx),numel(AP_idx));
end

% get mask of CPu + NAcc
function [striatum_mask, ccf_key] = get_str_mask(atlas_dir)
    % directories and atlas: see https://github.com/HoweLab/MultifiberLocalization
    ccf = load_tiffs_fast(fullfile(atlas_dir,'annotation_10_coronal.tif'));   
    ccf_key = load(fullfile(atlas_dir,'ccf_key.mat'),'labels','AP','ML','DV');
    
    % get striatum voxels
    is_str = ismember(ccf_key.labels.name,{'Nucleus accumbens','Caudoputamen'});
    str_values = unique(ccf_key.labels.value(is_str));

    % make striatum mask    
    striatum_mask = ismember(ccf,str_values);

    % save
    %bfsave(striatum_mask,'all_striatum_mask_10_coronal.tif','dimensionOrder','XYTZC','BigTiff',true,'Compression', 'LZW');
end
% calculate a bootstrap null distribution of local moran's I
% 1. shuffle data across recording locations
% 2. interpolate and smooth
% 3. calculate local moran's I
%
% Mai-Anh Vu, 2025
% updated Mai-Anh Vu, 2026/07/15

function output = local_morans_I_bootstrap_null(fiber_table,value_array,voxel_size,varargin)

%%%  parse optional inputs %%%
ip = inputParser;
ip.addParameter('n_it',10000);
ip.addParameter('AP_range',[]);
ip.addParameter('DV_range',[]);
ip.addParameter('ML_range',[]);
ip.addParameter('DV_sign',-1);
ip.addParameter('weight_matrix',ones(3,3,3));
ip.addParameter('return_null_data',0);
ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end

% initialize
output = struct;

% set up AP, ML, DV
if mode(sign(fiber_table.fiber_bottom_DV)) ~= DV_sign
    fiber_table.fiber_bottom_DV = -fiber_table.fiber_bottom_DV;
end
if isempty(AP_range)
    AP_range = [min(fiber_table.fiber_bottom_AP) max(fiber_table.fiber_bottom_AP)];
end
if isempty(ML_range)
    ML_range = [min(fiber_table.fiber_bottom_ML) max(fiber_table.fiber_bottom_ML)];
end
if isempty(DV_range)
    DV_range = [min(fiber_table.fiber_bottom_DV) max(fiber_table.fiber_bottom_DV)];
end

% fiber table formatting
if ischar(fiber_table.mouse)
    fiber_table.mouse = cellstr(fiber_table.mouse);
end
mice = unique(fiber_table.mouse);

% get interpolated map, interpolants, and original moran
tmp = get_activity_map_interp(value_array,fiber_table,voxel_size,...
    'AP_range',AP_range,'ML_range',ML_range,'DV_range',DV_range);
interp_struct = tmp.vol_01.interp_F;
output.moran = local_morans_I(tmp.vol_01.interp,'weight_matrix',weight_matrix);
clear tmp; 

% initialize these results
null_moran = nan(size(output.moran,1),size(output.moran,2),size(output.moran,3),n_it);


% shuffle non-nan values across recording locations within mouse
idx_it = nan(size(fiber_table,1),n_it);
for m = 1:numel(mice)
    mouse_idx = find(~isnan(value_array) & ismember(fiber_table.mouse,mice{m}));    
    for i = 1:n_it
        idx_it(mouse_idx,i) = vec(mouse_idx(randperm(numel(mouse_idx))));        
    end
end

% the loop
parfor i = 1:n_it
    % 1. shuffle
    these_vals = value_array(idx_it(:,i)); 
    % 2. interpolate
    tmp_sm = get_activity_map_interp(these_vals,fiber_table,voxel_size,...
        'AP_range',AP_range,'ML_range',ML_range,'DV_range',DV_range,...
        'interp_struct',interp_struct);
    tmp_sm = tmp_sm.vol_01.interp; % replace variable to help with memory
    % 3. local moran's I
    null_moran(:,:,:,i) = local_morans_I(tmp_sm,'weight_matrix',weight_matrix);
end

% output null moran data
if return_null_data == 1
    output.null_moran = null_moran;
end

% output mean and various prctiles
output.null_stats.mean = nanmean(int_results,4);
output.null_stats.prctile0p5 = prctile(int_results,0.5,4);
output.null_stats.prctile1 = prctile(int_results,1,4);
output.null_stats.prctile2p5 = prctile(int_results,2.5,4);
output.null_stats.prctile5 = prctile(int_results,5,4);
output.null_stats.prctile95 = prctile(int_results,95,4);
output.null_stats.prctile97p5 = prctile(int_results,97.5,4);
output.null_stats.prctile99 = prctile(int_results,99,4);
output.null_stats.prctile99p5 = prctile(int_results,99.5,4);

end



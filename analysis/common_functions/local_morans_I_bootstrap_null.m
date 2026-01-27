% calculate a bootstrap null distribution of local moran's I
% 1. shuffle data across recording locations
% 2. interpolate and smooth
% 3. calculate local moran's I
%
% Mai-Anh Vu, 2025
function output = local_morans_I_bootstrap_null(fib_tbl,vals,voxel_size,varargin)

%%%  parse optional inputs %%%
ip = inputParser;
ip.addParameter('n_it',10000);
ip.addParameter('AP_range',[]);
ip.addParameter('DV_range',[]);
ip.addParameter('ML_range',[]);
ip.addParameter('gaussian_smooth_sigma',2);
ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end

if isempty(AP_range)
    AP_range = prctile(fib_table.fiber_bottom_AP,[0 100]);
end
if isempty(ML_range)
    ML_range = prctile(fib_table.fiber_bottom_ML,[0 100]);
end
if isempty(DV_range)
    DV_range = prctile(fib_table.fiber_bottom_DV,[0 100]);
end

% get smooth map just for size of output
tmp = smooth_activity_map_interp_smoothed(vals,fib_tbl,voxel_size,...
    'AP_range',AP_range,'ML_range',ML_range,'DV_range',DV_range,...
    'gaussian_sigma',gaussian_smooth_sigma);
% preallocate output
sm_results = nan(size(tmp.smooth,1),size(tmp.smooth,2),size(tmp.smooth,3),n_it);

% shuffled indices to shuffle values across recording locations
idx_it = nan(size(fib_tbl,1),n_it);
for i = 1:n_it
    idx_it(:,i) = vec(randperm(size(fib_tbl,1)));
end

% the loop
output = struct;
parfor i = 1:n_it
    % 1. shuffle
    these_vals = vals(idx_it(:,i)); 
    % 2. smooth
    tmp = smooth_activity_map_interp_smoothed(these_vals,fib_tbl,voxel_size,...
        'AP_range',AP_range,'ML_range',ML_range,'DV_range',DV_range,...
        'gaussian_sigma',2);
    tmp = tmp.vol_01.smooth; % replace variable to help with memory
    % 3. local moran's I
    sm_results(:,:,:,i) = local_morans_I(tmp,'weight_matrix',ones(21,21,21));
end

% output mean and various prctiles
output.mean = nanmean(sm_results,4);
output.prctile0p5 = prctile(sm_results,0.5,4);
output.prctile1 = prctile(sm_results,1,4);
output.prctile2p5 = prctile(sm_results,2.5,4);
output.prctile5 = prctile(sm_results,5,4);
output.prctile95 = prctile(sm_results,95,4);
output.prctile97p5 = prctile(sm_results,97.5,4);
output.prctile99 = prctile(sm_results,99,4);
output.prctile99p5 = prctile(sm_results,99.5,4);

end



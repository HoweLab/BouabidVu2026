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
ip.addParameter('prctiles',[]);                 % list of percentiles (e.g., [95 97.5 99.5]) to evaluate; leave blank if none
ip.addParameter('it_sum',0);                    % return sum of each iteration? useful for moran neighborhood
ip.addParameter('batch_size',1000);             % save in chunks and read back in for memory saving; leave blank otherwise
ip.addParameter('save_null',[]);                % output path to save null distribution; leave blank if not saving

ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end

% initialize
output = struct;

% redundant, but seemingly  necessary for the parfor?
weight_matrix = weight_matrix;
n_it = n_it;
if isempty(batch_size)
    batch_size = n_it;
else
    batch_size = min([batch_size n_it]);
end

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

% shuffle non-nan values across recording locations within mouse
idx_it = nan(size(fiber_table,1),n_it);
for m = 1:numel(mice)
    mouse_idx = find(~isnan(value_array) & ismember(fiber_table.mouse,mice{m}));    
    for i = 1:n_it
        idx_it(mouse_idx,i) = vec(mouse_idx(randperm(numel(mouse_idx))));        
    end
end
output.null_shuffle = idx_it;

% if we're saving out the null distribution, start the file
if ~isempty(save_null)
    output_null = matfile(save_null, 'Writable', true);
    output_null.null_moran = zeros(size(output.moran,1),size(output.moran,2),size(output.moran,3),n_it,'single'); % single for memory
end

% if we're keeping prctiles, keep an upper buffer of max values
if ~isempty(prctiles)
    n_to_keep =  ceil((100-min(prctiles))/100*n_it);
    upper_buffer = [];
end

% preallocate if we're taking the sum of each iteration iteration
if it_sum == 1
    output.null_stats.it_sum = nan(n_it,1);
end

% the loop
disp('     starting null')
n_chunk = ceil(n_it/batch_size);
for j = 1:n_chunk
    null_moran = nan(size(output.moran,1),size(output.moran,2),size(output.moran,3),batch_size);
    parfor i = 1:batch_size
        % 1. shuffle
        these_vals = value_array(idx_it(:,(j-1)*batch_size+i)); 
        % 2. interpolate
        tmp_sm = get_activity_map_interp(these_vals,fiber_table,voxel_size,...
            'AP_range',AP_range,'ML_range',ML_range,'DV_range',DV_range,...
            'interp_struct',interp_struct);
        tmp_sm = tmp_sm.vol_01.interp; % replace variable to help with memory
        % 3. local moran's I
        null_moran(:,:,:,i) = local_morans_I(tmp_sm,'weight_matrix',weight_matrix);
    end
    
    % upper buffer
    if ~isempty(prctiles)
        upper_buffer = cat(4,upper_buffer,null_moran);
        upper_buffer = sort(upper_buffer,4,'descend');
        upper_buffer = upper_buffer(:,:,:,1:n_to_keep);
    end
    
    % iteration sum
    if it_sum == 1
        output.null_stats.it_sum((j-1)*batch_size+(1:batch_size)) = ...
            vec(nansum(null_moran,[1 2 3]));
    end
    
    % write
    output_null.null_moran(:,:,:,(j-1)*batch_size+(1:batch_size)) = single(null_moran); % cast to single   
    disp(['     done: ' num2str(batch_size*j)])
end

% output prctiles
if ~isempty(prctiles)
    upper_buffer = sort(upper_buffer,4,'ascend');
    for p = 1:numel(prctiles)
        p_field = (['prctile' strrep(sprintf('%2.02f',prctiles(p)),'.','p')]);
        p_idx = round(prctiles(p)/100*n_it-n_it+n_to_keep+1);
        output.null_stats.(p_field) = upper_buffer(:,:,:,p_idx);
    end
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_moran_neighborhood:
%
% Try a handful of radii and a quick null to get an empirical best
% neighborhood size

function moran_neighborhood = get_moran_neighborhood_width(fiber_table,value_array,voxel_size,widths,varargin)

%%%  parse optional inputs %%%
ip = inputParser;
ip.addParameter('n_it',500);
ip.addParameter('AP_range',[]);
ip.addParameter('DV_range',[]);
ip.addParameter('ML_range',[]);
ip.addParameter('DV_sign',-1);
ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end

% initialize
moran_neighborhood = struct;
moran_neighborhood.width = vec(widths);
moran_neighborhood.width_z = nan(numel(widths),1);


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

for w = 1:numel(widths)
    this_width = widths(w);
    disp(['     neighborhood width: ' num2str(this_width)])
    weight_matrix = ones(this_width,this_width,this_width);
    tmp = local_morans_I_bootstrap_null(fiber_table,value_array,voxel_size,...
        'AP_range',AP_range,'ML_range',ML_range,'DV_range',DV_range,'n_it',n_it,...
        'weight_matrix',weight_matrix,'it_sum',1);
    moran_neighborhood.width_z(w) = (nansum(tmp.moran(:)) - mean(tmp.null_stats.it_sum))/std(tmp.null_stats.it_sum); 
end


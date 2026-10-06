function output = get_contour_volumes(input_volume,varargin)

%%%  parse optional inputs %%%
ip = inputParser;
ip.addParameter('contour_thresh',[]); 
ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end

% contour thresh
if isempty(contour_thresh)
    contour_thresh = prctile(input_volume(:),[0:10:100]); % use 0.01 instead of 0 for the first one
end

% initialize
output = struct;
output.contour.thresh = contour_thresh;

% max
[max_val,max_idx] = nanmax(input_volume(:));       
[i,j,k] = ind2sub(size(input_volume),max_idx);    
output.max.val = max_val;
output.max.idx = max_idx;
output.max.row_col_slice = [i j k];
     
% min
[min_val,min_idx] = nanmin(input_volume(:));
[i,j,k] = ind2sub(size(input_volume),min_idx);    
output.min.val = min_val;
output.min.idx = min_idx;
output.min.row_col_slice = [i j k];


% 3D contour volumes

% simple discretize
output.contour.discretize = discretize(input_volume,contour_thresh);
replace_vals = movmean(contour_thresh,2);
replace_vals = replace_vals(2:end);
output.contour.discretize(~isnan(output.contour.discretize)) = ...
    replace_vals(output.contour.discretize(~isnan(output.contour.discretize)));

% less than
all_thresh_vol = zeros(size(input_volume));
for t = 2:numel(contour_thresh)
    
    % threshold
    this_thresh_vol = input_volume <= contour_thresh(t);
    
    % keep the biggest one that has the minimum in it
    conn_vol = bwconncomp(this_thresh_vol);
    conn_vol = conn_vol.PixelIdxList;
    conn_vol = conn_vol(cellfun(@(x) ismember(output.min.idx,x),conn_vol));    
    [~,biggest_vol]=max(cellfun(@numel,conn_vol));

    % contour volume
    this_thresh_vol = zeros(size(this_thresh_vol));    
    this_thresh_vol(conn_vol{biggest_vol}) = 1;
    
    % add
    all_thresh_vol = all_thresh_vol + this_thresh_vol;
end
replace_vals = fliplr(contour_thresh(2:end));
all_thresh_vol_repl = nan(size(all_thresh_vol));
all_thresh_vol_repl(all_thresh_vol>0) = replace_vals(all_thresh_vol(all_thresh_vol>0));
output.contour.less_than = all_thresh_vol_repl;

% greater than
all_thresh_vol = zeros(size(input_volume));
for t = (numel(contour_thresh)-1):-1:1
    
    % threshold
    this_thresh_vol = input_volume >= contour_thresh(t);
    
    % keep the biggest one that has the maximum in it
    conn_vol = bwconncomp(this_thresh_vol);
    conn_vol = conn_vol.PixelIdxList;
    conn_vol = conn_vol(cellfun(@(x) ismember(output.max.idx,x),conn_vol));
    [~,biggest_vol]=max(cellfun(@numel,conn_vol));

    % contour volume
    this_thresh_vol = zeros(size(this_thresh_vol));    
    this_thresh_vol(conn_vol{biggest_vol}) = 1;
    
    % add
    all_thresh_vol = all_thresh_vol + this_thresh_vol;
end
replace_vals = contour_thresh(1:end-1);
all_thresh_vol_repl = nan(size(all_thresh_vol));
all_thresh_vol_repl(all_thresh_vol>0) = replace_vals(all_thresh_vol(all_thresh_vol>0));
output.contour.greater_than = all_thresh_vol_repl;




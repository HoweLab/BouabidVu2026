% function M = local_morans_I(data,varargin)
%
% input: data = a 3D matrix of values (nans will be handled appropriately)
%   - fiber data can be transformed to this format via interpolation or
%   interpolation plus smoothing
%
% optional input: weight_matrix, determines the neighborhood over which
% local Moran's I will be calculated
%
% Calculate local morans I for 3D gridded data based on a neighbor weight
% matrix. If none is given, the default is the 3d equivalent of "queen
% weights" of a single layer
%
% see https://storymaps.arcgis.com/stories/5b26f25bb81a437b89003423505e2f71
%
% VECTORIZED VERSION: mathematically equivalent to the original
% nested-loop implementation, but replaces the per-voxel loop + sub2ind
% + get_m() calls with a handful of N-D convolutions (convn). This is
% pure base-MATLAB (no toolbox dependency beyond nanmean, which the
% original already required).
%
% NOTE ON EQUIVALENCE: the original loop performs *correlation*
% (neighbor offsets applied directly), while convn performs true
% *convolution* (kernel flipped in every dimension). This version flips
% the kernel before calling convn to restore correlation semantics, so
% results match the original exactly -- including for asymmetric,
% user-supplied weight_matrix inputs, not just the symmetric default.
%
% original by Mai-Anh Vu, 2025/09/11, see bottom.
% vectorized rewrite by Claude Sonnet5 2026/07/14
% validated and compared by Mai-Anh, 2026/07/15
function M = local_morans_I(data,varargin)

%%% parse optional inputs %%%
ip = inputParser;
ip.addParameter('weight_matrix',ones(3,3,3));
ip.parse(varargin{:});
weight_matrix = ip.Results.weight_matrix;

% zero out the center (self-weight), same as original
W = weight_matrix;
W(ceil(numel(W)/2)) = 0;

% if the data are embedded shells of nans, trim to the bounding box for speed
row_keep = [find(sum(~isnan(data),[2 3]) > 0,1,'first'),...
    find(sum(~isnan(data),[2 3]) > 0,1,'last')];
col_keep = [find(sum(~isnan(data),[1 3]) > 0,1,'first'),...
    find(sum(~isnan(data),[1 3]) > 0,1,'last')];
slice_keep = [find(sum(~isnan(data),[1 2]) > 0,1,'first'),...
    find(sum(~isnan(data),[1 2]) > 0,1,'last')];

M = nan(size(data));

sub = data(row_keep(1):row_keep(2),col_keep(1):col_keep(2),slice_keep(1):slice_keep(2));

% mask of valid (non-nan) voxels, and mean-centered data with nans -> 0
% (nans contribute exactly 0 to any weighted sum below, matching the
% original's nansum-based masking behavior)
mask = double(~isnan(sub));
z = sub - nanmean(sub(:));
z(isnan(z)) = 0;

% correlation kernels (flip convn's kernel to undo its internal flip)
flipdims = @(A) A(end:-1:1,end:-1:1,end:-1:1);
Wcorr    = flipdims(W);           % actual (possibly non-binary) weights
Wbin     = double(W ~= 0);        % binary indicator of "counts as a neighbor"
Wbincorr = flipdims(Wbin);

% --- numerator / denominator for the weighted neighbor average ---
% (equivalent to nansum(mean_diff_sub .* normalized_weight_matrix) per voxel)
num = convn(z,    Wcorr,    'same');
den = convn(mask, Wcorr,    'same');
weighted_avg = num ./ den;   % NaN/Inf where a voxel has zero valid weighted neighbors (edge case, matches original's degenerate behavior)

% --- Si2: local variance term over valid, kernel-nonzero neighbors ---
sq_num = convn(z.^2, Wbincorr, 'same');
cnt    = convn(mask,  Wbincorr, 'same');
Si2 = sq_num ./ (cnt - 1);

% --- assemble local Moran's I ---
lm = z .* weighted_avg ./ Si2;

M(row_keep(1):row_keep(2),col_keep(1):col_keep(2),slice_keep(1):slice_keep(2)) = lm;
M(isnan(data)) = nan; % put nans back
end


% % ORIGINAL VERSION BY MAI-ANH, basis for Claude vectorized rewrite above
% % function M = local_morans_I(data,varargin) 
% % 
% % input: data = a 3D matrix of values (nans will be handled appropriately)
% %   - fiber data can be transformed to this format via interpolation or
% %   interpolation plus smoothing
% %
% % optional input: weight_matrix, determines the neighborhood over which
% % local Moran's I wil be calculated
% %
% % Calculate local morans I for 3D gridded data based on a neighbor weight
% % matrix. If none is given, the default is the 3d equivalent of "queen
% % weights" of a single layer
% %
% % see https://storymaps.arcgis.com/stories/5b26f25bb81a437b89003423505e2f71
% %
% %
% % Mai-Anh Vu, 2025/09/11
% 
% function M = local_morans_I(data,varargin) 
% 
% %%%  parse optional inputs %%%
% ip = inputParser;
% ip.addParameter('weight_matrix',ones(3,3,3));
% ip.parse(varargin{:});
% for j=fields(ip.Results)'
%     eval([j{1} '=ip.Results.' j{1} ';']);
% end
% 
% % if the data are embedded shells of nans, let's do this for speed
% row_keep = [find(sum(~isnan(data),[2 3]) > 0,1,'first'),...
%     find(sum(~isnan(data),[2 3]) > 0,1,'last')];
% col_keep = [find(sum(~isnan(data),[1 3]) > 0,1,'first'),...
%     find(sum(~isnan(data),[1 3]) > 0,1,'last')];
% slice_keep = [find(sum(~isnan(data),[1 2]) > 0,1,'first'),...
%     find(sum(~isnan(data),[1 2]) > 0,1,'last')];
% M = nan(size(data));
% 
% % now data is just the non-nan cube within data
% data = data(row_keep(1):row_keep(2),col_keep(1):col_keep(2),slice_keep(1):slice_keep(2));
% 
% % difference of all data from the mean
% data_minus_mean = data-nanmean(data(:));
% 
% % weight matrix
% W = vec(weight_matrix); % vector for convenience
% W(ceil(numel(W)/2)) = 0; % zero out middle
% % note: it will get normalized in the sub function
% 
% % mini sliding cubes of relative indices (vectorized into single col)
% cube_row = repmat(vec(1:size(weight_matrix,1)),...
%     1,size(weight_matrix,2),size(weight_matrix,3));
% cube_row = cube_row-median(cube_row);
% cube_col = vec(permute(cube_row,[2 1 3]));
% cube_slice = vec(permute(cube_row,[3 2 1]));
% cube_row = vec(cube_row);
% 
% % loop over voxels
% lm = nan(size(data));
% for i = 1:size(data,1)
%     this_cube_row = cube_row+i;
%     this_cube_row(this_cube_row<1 | this_cube_row>size(data,1)) = nan;
%     for j = 1:size(data,2)
%         this_cube_col = cube_col+j;
%         this_cube_col(this_cube_col<1 | this_cube_col>size(data,2)) = nan;
%         for k = 1:size(data,3)
%             this_cube_slice = cube_slice+k;
%             this_cube_slice(this_cube_slice<1 | this_cube_slice>size(data,3)) = nan;
%             this_ind = sub2ind(size(data),this_cube_row,this_cube_col,this_cube_slice);
%             this_cube = nan(size(this_ind));
%             this_cube(~isnan(this_ind)) = data_minus_mean(this_ind(~isnan(this_ind)));
%             lm(i,j,k) = get_m(this_cube,W);
%         end
%     end
% end
% 
% % put back into bigger output
% M(row_keep(1):row_keep(2),col_keep(1):col_keep(2),slice_keep(1):slice_keep(2)) = lm;
% 
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%% the actual calculation
% function local_m = get_m(mean_diff_sub,weight_matrix)
%     
%     % for convenience. might already be in this format, but just in case not
%     mean_diff_sub = vec(mean_diff_sub); 
%     weight_matrix(isnan(mean_diff_sub)) = 0; % handle nans
%     weight_matrix = vec(weight_matrix)/sum(vec(weight_matrix));
%     weight_matrix_ones = double(weight_matrix>0); 
%     
%     % now the the calculation
%     xi_minus_X = mean_diff_sub(ceil(numel(mean_diff_sub)/2)); % center value
%     Si2 = nansum((mean_diff_sub.*weight_matrix_ones).^2)/...
%         (sum(weight_matrix_ones)-1);
%     local_m = xi_minus_X/Si2*nansum(mean_diff_sub.*weight_matrix);
%     
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%% quick function to turn matrices into column vector
% function y = vec(x)
%     y = x(:);
% end
% 


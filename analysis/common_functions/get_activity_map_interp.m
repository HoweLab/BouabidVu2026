% Mai-Anh Vu, 2025/09/10
% notes on inputs
%   value_array has 1 value for every entry in the fiber_table
%   voxel_size is in mm
%
% edited 26/07/15 to handle pseudo-replication, and remove smoothing
% edited 26/10/08 to handle distance weighting and boundary extrapolation
% (needs MATLAB 2024a or later)


function output = get_activity_map_interp(value_array,fiber_table,voxel_size,varargin)

%%%  parse optional inputs %%%
ip = inputParser;
ip.addParameter('AP_range',[]);
ip.addParameter('ML_range',[]);
ip.addParameter('DV_range',[]);
ip.addParameter('DV_sign',-1); 
ip.addParameter('interp_struct',[]);
ip.addParameter('incl_plot_info',false);
ip.addParameter('incl_projections',false);
ip.addParameter('dist_weighting',true);
ip.addParameter('mouse_dist_sigma',0.15);
ip.addParameter('hard_dist_cutoff',[]);
ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
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

% set up the meshgrid space
x = ML_range(1):voxel_size:ML_range(2);
y = DV_range(1):voxel_size:DV_range(2);
z = AP_range(1):voxel_size:AP_range(2);
[xx,yy,zz] = meshgrid(x,y,z);

% all coordinates
fib_coord = [fiber_table.fiber_bottom_ML,...
    fiber_table.fiber_bottom_DV,...
    fiber_table.fiber_bottom_AP];

% other setup
if isempty(hard_dist_cutoff)
    hard_dist_cutoff = 3*mouse_dist_sigma;
end

% make sure mouse name is cellstr
if ischar(fiber_table.mouse)
    fiber_table.mouse = cellstr(fiber_table,mouse);
end

% initialize
output = struct;
output.info.dist_weighting = dist_weighting; % if false, support will be n_mice

% loop over columns
for i = 1:size(value_array,2)
        
    % now loop over mice    
    vals = value_array(:,i);        % these vals
    idx = ~isnan(vals);             % get non-nan idx
    mice = unique(fiber_table.mouse(idx));  % mice for this value array        
    
    % initialize 4D matrices to track mouse interp maps, support weights, etc
    mouse_interp = nan(size(xx,1),size(xx,2),size(xx,3),numel(mice));
    mouse_weight = zeros(size(xx,1),size(xx,2),size(xx,3),numel(mice));
    dnear = inf(size(xx,1),size(xx,2),size(xx,3),numel(mice));
    
    for m = 1:numel(mice)
        mouse_idx = ~isnan(vals) & ismember(fiber_table.mouse,mice(m));
        if sum(mouse_idx)>=4
        
            % natural neighbor interpolation        
            v = version;
            if isempty(interp_struct)
                if str2double(v(1:2)) < 24
                    F = scatteredInterpolant(fib_coord(mouse_idx,:),vals(mouse_idx),'natural','nearest');
                else
                    F = scatteredInterpolant(fib_coord(mouse_idx,:),vals(mouse_idx),'natural','boundary');
                end
            else % if we already have the interpolant, just update the Values
                F = interp_struct.(get_mouse_field(mice{m}));
                F.Values = vals(mouse_idx);
            end
    
            % interpolated volume
            mouse_interp(:,:,:,m) = F(xx,yy,zz);
            output.(['vol_' sprintf('%02d',i)]).interp_F.(get_mouse_field(mice{m})) = F; % mouse interpolant
          
            % distance weighting
            if dist_weighting
                [~, d] = knnsearch(fib_coord, [xx(:) yy(:) zz(:)]);   % distance to this mouse's nearest fiber
                dnear(:,:,:,m) = reshape(d,size(dnear,1),size(dnear,2),size(dnear,3));
                w = exp(-0.5 * (d / mouse_dist_sigma).^2);
                w(d > hard_dist_cutoff) = 0;
                w(isnan(mouse_interp(:,:,:,m))) = 0;
                mouse_weight(:,:,:,m) = reshape(w,size(mouse_weight,1),size(mouse_weight,2),size(mouse_weight,3));
            end
        end
    end    
    
    % group mean & support for each voxel
    if ~dist_weighting
         output.(['vol_' sprintf('%02d',i)]).interp = reshape(mean(all_vol,2,'omitnan'),numel(y),numel(x),numel(z));
         output.(['vol_' sprintf('%02d',i)]).support = reshape(sum(~isnan(mouse_interp),2),numel(y),numel(x),numel(z));
    else
        group_weight = sum(mouse_weight,4);
        mouse_interp_0 = mouse_interp;
        mouse_interp_0(mouse_weight==0) = 0; % avoid nan*0 = nan
        group_vol = sum(mouse_weight.*mouse_interp_0, 2)./group_weight;
        group_vol(group_weight==0) = nan; % nan        
        output.(['vol_' sprintf('%02d',i)]).interp = group_vol;
        output.(['vol_' sprintf('%02d',i)]).support = group_weight;
    end
    
    if incl_projections
        % projections for convenience
        output.(['vol_' sprintf('%02d',i)]).axial.mean_projection = permute(nanmean(output.(['vol_' sprintf('%02d',i)]).interp,1),[3 2 1]);
        output.(['vol_' sprintf('%02d',i)]).sagittal.mean_projection = permute(nanmean(output.(['vol_' sprintf('%02d',i)]).interp,2),[1 3 2]);
        output.(['vol_' sprintf('%02d',i)]).coronal.mean_projection = permute(nanmean(output.(['vol_' sprintf('%02d',i)]).interp,3),[1 2 3]);
    end
end

% some basic info
if incl_plot_info
    output.info.voxel_size = voxel_size;
    output.info.dimension_order = {'DV','ML','AP';'y','x','z'};
    output.info.ML = x;
    output.info.AP = z;
    output.info.DV = y;
    output.info.all_ML = xx;
    output.info.all_AP = zz;
    output.info.all_DV = yy;

    output.info.axial.flatten_dim = 1;
    output.info.axial.permute = [3 2 1];
    output.info.axial.x = output.info.ML;
    output.info.axial.y = output.info.AP;
    output.info.axial.x_dir = 'normal';
    output.info.axial.y_dir = 'normal';

    output.info.sagittal.flatten_dim = 2;
    output.info.sagittal.permute = [1 3 2];
    output.info.sagittal.x = output.info.AP;
    output.info.sagittal.y = output.info.DV;
    output.info.sagittal.x_dir = 'reverse';
    output.info.sagittal.y_dir = 'normal';

    output.info.coronal.flatten_dim = 3;
    output.info.coronal.permute = [1 2 3];
    output.info.coronal.x = output.info.ML;
    output.info.coronal.y = output.info.DV;
    output.info.coronal.x_dir = 'normal';
    output.info.coronal.y_dir = 'normal';
end

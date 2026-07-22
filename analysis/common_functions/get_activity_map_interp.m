% Mai-Anh Vu, 2025/09/10
% notes on inputs
%   value_array has 1 value for every entry in the fiber_table
%   voxel_size is in mm
%
% edited 26/07/15 to handle pseudoreplication, and remove smoothing


function output = get_activity_map_interp(value_array,fiber_table,voxel_size,varargin)

%%%  parse optional inputs %%%
ip = inputParser;
ip.addParameter('AP_range',[]);
ip.addParameter('ML_range',[]);
ip.addParameter('DV_range',[]);
ip.addParameter('DV_sign',-1); 
ip.addParameter('interp_struct',[]);
ip.addParameter('incl_plot_info',0);
ip.addParameter('incl_projections',0);
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

% make sure mouse name is cellstr
if ischar(fiber_table.mouse)
    fiber_table.mouse = cellstr(fiber_table,mouse);
end

% loop over columns
for i = 1:size(value_array,2)
        
    % now loop over mice    
    vals = value_array(:,i);        % these vals
    idx = ~isnan(vals);             % get non-nan idx
    mice = unique(fiber_table.mouse(idx));  % mice for this value array
    
    % initialize 4D matrix to track mouse interp maps
    mouse_interp = nan(size(xx,1),size(xx,2),size(xx,3),numel(mice));
    
    for m = 1:numel(mice)
        mouse_idx = ~isnan(vals) & ismember(fiber_table.mouse,mice(m));
        
        % natural neighbor interpolation
        if isempty(interp_struct)
            F = scatteredInterpolant( ...
                fiber_table.fiber_bottom_ML(mouse_idx),...
                fiber_table.fiber_bottom_DV(mouse_idx),...
                fiber_table.fiber_bottom_AP(mouse_idx),...
                vals(mouse_idx),...
                'natural','none');
        else % if we already have the interpolant, just update the Values
            F = interp_struct.(mice{m});
            F.Values = vals(mouse_idx);
        end
        % interpolated volume
        mouse_interp(:,:,:,m) = F(xx,yy,zz);
        output.(['vol_' sprintf('%02d',i)]).interp_F.(mice{m}) = F;
    end    
       
    % group mean & #mice contributing to each voxel
    output.(['vol_' sprintf('%02d',i)]).interp = nanmean(mouse_interp,4);
    output.(['vol_' sprintf('%02d',i)]).n_mice = sum(~isnan(mouse_interp),4);
    
    if incl_projections == 1
        % projections for convenience
        output.(['vol_' sprintf('%02d',i)]).axial.mean_projection = permute(nanmean(output.(['vol_' sprintf('%02d',i)]).interp,1),[3 2 1]);
        output.(['vol_' sprintf('%02d',i)]).sagittal.mean_projection = permute(nanmean(output.(['vol_' sprintf('%02d',i)]).interp,2),[1 3 2]);
        output.(['vol_' sprintf('%02d',i)]).coronal.mean_projection = permute(nanmean(output.(['vol_' sprintf('%02d',i)]).interp,3),[1 2 3]);
    end
end

% some basic info
if incl_plot_info==1
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

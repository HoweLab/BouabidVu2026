% Mai-Anh Vu, 2025/09/10
% notes on inputs
%   value_array has 1 value for every entry in the fiber_table
%   voxel_size is in mm

function output = smooth_activity_map_interp_smoothed(value_array,fiber_table,voxel_size,varargin)

%%%  parse optional inputs %%%
ip = inputParser;
ip.addParameter('AP_range',[]);
ip.addParameter('ML_range',[]);
ip.addParameter('DV_range',[]);
ip.addParameter('DV_sign',-1); 
ip.addParameter('gaussian_sigma',1); 
ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end

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


x = ML_range(1):voxel_size:ML_range(2);
y = DV_range(1):voxel_size:DV_range(2);
z = AP_range(1):voxel_size:AP_range(2);
[xx,yy,zz] = meshgrid(x,y,z);

output.info.voxel_size = voxel_size;
output.info.dimension_order = {'DV','ML','AP';'y','x','z'};
output.info.ML = x;
output.info.AP = z;
output.info.DV = y;
output.info.all_ML = xx;
output.info.all_AP = zz;
output.info.all_DV = yy;


for i = 1:size(value_array,2)
    vals = value_array(:,i);    
    idx = ~isnan(vals);
    F = scatteredInterpolant( ...
        fiber_table.fiber_bottom_ML(idx),...
        fiber_table.fiber_bottom_DV(idx),...
        fiber_table.fiber_bottom_AP(idx),...
        vals(idx),...
        'natural','none');
    vv = F(xx,yy,zz);
    
    
    % see here: https://stackoverflow.com/questions/72722767/how-to-pad-an-irregularly-shaped-matrix-in-matlab
    % we need to pad this so that the gaussian filter doesn't bring the edges in
    
    % Track the "inner" indices, where values are defined
    inner = (~isnan(vv));
    B = vv;                 % Copy A so we don't change it
    for j = 1:2*ceil(2*gaussian_sigma)+1
        B(~inner) = 0;                  % Replace NaN with 0 so that convolutions work OK
        % First dilate the inner region by one element, taking the average of
        % neighbours which are up/down/left/right (no diagonals). This is required
        % to avoid including interior points (which only touch diagonally) in the
        % averaging. These can be considered the "cardinal neighbours"
        kernel = cat(3,cat(3,[0 0 0; 0 1 0; 0 0 0],...
            [0 1 0 ; 1 0 1; 0 1 0]),[0 0 0; 0 1 0; 0 0 0]); % Cardinal directions in 3x3 stencil
        s = convn(B,kernel,'same');      % 2D convolution to get sum of neighbours
        n = convn(inner,kernel,'same');  % 2D convolution to get count of neighbours
        s(inner) = 0;                    % Zero out the inner region
        s = s./n;                        % Get the mean of neighbours

        % Second, dilate the inner region but including the mean from all
        % directions. This lets us handle convex corners in the image
        s2 = convn(B,ones(3,3,3),'same');     % Sum of neighbours (and self, doesn't matter)
        n = convn(inner,ones(3,3,3),'same');  % Count of neighbours (self=0 for dilated elems)
        s2 = s2./n;                       % Get the mean of neighbours

        % Finally piece together the 3 matrices:
        out = s2;                       % Start with outmost dilation inc. corners
        out(~isnan(s)) = s(~isnan(s));  % Override with inner dilation for cardinal neighbours
        out(inner) = B(inner);          % Override with original inner data
        
        % update
        B = out;
        inner = ~isnan(B);
    end

    output.(['vol_' sprintf('%02d',i)]).interp = vv;
    output.(['vol_' sprintf('%02d',i)]).intermed = B;
    output.(['vol_' sprintf('%02d',i)]).smooth = imgaussfilt3(B,gaussian_sigma,'FilterDomain','spatial');    
    
    % projections for convenience
    output.(['vol_' sprintf('%02d',i)]).axial.mean_projection = permute(nanmean(output.(['vol_' sprintf('%02d',i)]).smooth,1),[3 2 1]);
    output.(['vol_' sprintf('%02d',i)]).sagittal.mean_projection = permute(nanmean(output.(['vol_' sprintf('%02d',i)]).smooth,2),[1 3 2]);
    output.(['vol_' sprintf('%02d',i)]).coronal.mean_projection = permute(nanmean(output.(['vol_' sprintf('%02d',i)]).smooth,3),[1 2 3]);
    
    % plot_info.axial.flatten_dim = 1;
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



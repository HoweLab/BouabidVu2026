% function output = get_volume_projection(this_map,orientation,varargin)
%
% optional inputs
% mask: mask volume; values outside mask will be assigned nan
% proj: what projection operation; default = 'mean'; other options = 'max','min'
%
% Mai-Anh Vu, 2025
function output = get_volume_projection(this_map,orientation,varargin)

%%%  parse optional inputs %%%
ip = inputParser;
% directories and atlas: see https://github.com/HoweLab/MultifiberLocalization    
ip.addParameter('mask',[]);
ip.addParameter('mask_replace',0);
ip.addParameter('proj','mean'); % projection operation. other options = max, min
ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end

output = this_map;

% if we're masking
if ~isempty(mask)
    output(mask==0) = mask_replace;
end
 
% permute orientation
permutes.axial = [3 2 1];
permutes.sagittal = [1 3 2];
permutes.coronal = [1 2 3];

% dimension for operation
op_dim.axial = 1;
op_dim.sagittal = 2;
op_dim.coronal = 3;

if strcmp(proj,'mean')
    output = permute(nanmean(output,op_dim.(orientation)),permutes.(orientation));
elseif strcmp(proj,'max')
    output = permute(nanmax(output,[],op_dim.(orientation)),permutes.(orientation));
elseif strcmp(proj,'min')
    output = permute(nanmin(output,[],op_dim.(orientation)),permutes.(orientation));
else
    output = [];
end

    
    

% function get_dropoff_from_extrema(input_volume,peak_def,str,map_info,varargin)
% 
% This function calculates the change in value as a function of distance
% from peak, defined as either the max or min of the map
%
% Inputs:
% input_volume: interpolated volume
% peak_def: 'max' or 'min'

function output = get_dropoff_from_extrema(input_volume,peak_def,str,map_info,varargin)

% mask
input_volume(str.striatum_mask==0) = nan;

% get peak info
if strcmp(peak_def,'min')
    % min
    [peak_val,peak_idx] = nanmin(input_volume(:));
elseif strcmp(peak_def,'max')
    % max
    [peak_val,peak_idx] = nanmax(input_volume(:));       
end

% peak coordinate (mm from Bregma)
[i,j,k] = ind2sub(size(input_volume),peak_idx);        
peak_ML = map_info.ML(j);
peak_DV = map_info.DV(i);
peak_AP = map_info.AP(k);

% now get all voxels
all_vals = input_volume(:)-peak_val;
all_ML = abs(map_info.all_ML(:)-peak_ML);
all_DV = abs(map_info.all_DV(:)-peak_DV);
all_AP = abs(map_info.all_AP(:)-peak_AP);

% get rid of nans
all_ML(isnan(all_vals)) = [];
all_DV(isnan(all_vals)) = [];
all_AP(isnan(all_vals)) = [];
all_vals(isnan(all_vals)) = [];
% 
% figure
% subplot(1,3,1)
% plot(all_AP,all_vals,'.k')
% title('AP')
% subplot(1,3,2)
% plot(all_ML,all_vals,'.k')
% title('ML')
% subplot(1,3,3)
% plot(all_DV,all_vals,'.k')
% title('DV')

% model
tbl = table(all_vals,all_AP,all_ML,all_DV,'VariableNames',{'val','AP','ML','DV'});
mdl_spec = 'val ~ AP + ML + DV';
mdl = fitglm(tbl,mdl_spec);

% output
output = struct;
output.peak = peak_val;
output.peak_index = peak_idx;
output.peak_AP = peak_AP;
output.peak_ML = peak_ML;
output.peak_DV = peak_DV;
output.mdl = mdl;

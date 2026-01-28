% function output = raster_transient_correlation(ref_raster,corr_raster,polarity,varargin)
%
%
% Example use case: for unpredicted reward: at what temporal latency is 
% trial-by-trial ACh correlated with trial-by-trial DA peak magnitude
%
%
% inputs: 
% ref_raster (mxnxp, m=timepoints, n=ROIs, p=trials): raster in which I 
%       care about the peak or trough (e.g. DA in the above example)
% corr_raster (mxn, m=timepoints, n=ROIs, p=trials): raster in which I'll 
%       evaluate the correlation at different temporal lags 
% input_idx: the indices aligning with the rows of the rasters, i.e., -9:27
% polarty: 1 (peak) or -1 (trough)
%
% optional inputs:
% peak_idx_of_int: which indices of the raster to consider for the
%       transient (default = all)
% corr_idx_of_int: which indices of the raster to consider for the
%       correlation (default = all)
% local_transient_window: the window around the max (or min) of the
%       triggered average over which to consider the trial_wise max (or
%       min). it is a symmetric window; i.e., 18 = +/-9 
% center_transient: if set to 1, each trial will center the peak or trough
%       of the transient as t=0, otherwise the average peak or trough time
%       will be used as t=0
%
% output is a struct with lots of information (see bottom for details)
% 
% 
%
% Mai-Anh Vu
% 2025/09/16



function output = raster_transient_correlation(ref_raster,corr_raster,polarity,varargin)

%%%  parse optional inputs %%%
ip = inputParser;
% the indices aligning with the rows of the rasters, i.e., -9:27; default
% 1:_
ip.addParameter('input_idx',[]); 
% which indices of the raster to consider for the transient
ip.addParameter('ref_idx_of_int',[]); 
% which indices of the raster to consider for the correlation
ip.addParameter('corr_idx_of_int',1); 
% window around max (or min) of triggered average over which to consider 
% trial-wise max (or min). a window of 1 consider +/-1 index around the
% triggered average max index
ip.addParameter('local_transient_window',[]);
% default: 1 = center the trial-wise transient as t=0; else 0 = center the
% triggered-average transient
ip.addParameter('center_trial_transient',0); 
ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end

% indices of interest
if isempty(input_idx)
    input_idx = 1:size(ref_raster,1);
end
if isempty(ref_idx_of_int)
    ref_idx_of_int = 1:size(ref_raster,1);
end
if isempty(corr_idx_of_int)
    corr_idx_of_int = 1:size(corr_raster,1);
end


% nan out the ref_raster indices we don't care about
ref_raster_orig = ref_raster;
ref_raster(~ismember(input_idx,ref_idx_of_int),:,:) = nan;
 

% get triggered average max/min to anchor timing of finding local max/min
if polarity == 1
    [~,avg_tr_idx] = nanmax(nanmean(ref_raster,3),[],1); 
else
    [~,avg_tr_idx] = nanmin(nanmean(ref_raster,3),[],1); 
end

% get local max nearest triggered average max
%
% here's what's happening here: i'm convolving the trg avg max with 
% something like a gaussian... except it is a linear peak thing, just for 
% simplicity. that way, i can element-wise multiply that with the local 
% maxes/mins (ones) and find the max for a built-in tiebreaker. Then I'll
% add the DFF values (they're all <1, i.e., too small to change the 
% outcome in clear winners, just if 2 local maxes are the same distance 
% from the trg avg max/min)
if isempty(local_transient_window)
    local_transient_window = floor(sum(~isnan(ref_raster(1,:,1)))/2);
end
% make a new raster centered on the triggered average max/min
loc_idx_i = repmat(vec(-local_transient_window:local_transient_window),...
    1,size(ref_raster,2),size(ref_raster,3)) + ...
    repmat(avg_tr_idx,local_transient_window*2+1,1,size(ref_raster,3));
loc_idx_i(loc_idx_i<1 | loc_idx_i>size(ref_raster,1)) = nan;
loc_idx_j = repmat(1:size(ref_raster,2),local_transient_window*2+1,1,size(ref_raster,3));
loc_idx_k = repmat(permute(1:size(ref_raster,3),[1 3 2]),local_transient_window*2+1,size(ref_raster,2));
loc_idx_ind = sub2ind(size(ref_raster), loc_idx_i,loc_idx_j,loc_idx_k);
ref_raster_local = nan(size(loc_idx_ind));
ref_raster_local(~isnan(loc_idx_ind)) = ref_raster_orig(loc_idx_ind(~isnan(loc_idx_ind)));
if polarity == 1
    local_tr = ref_raster_local>0 & islocalmax(ref_raster_local,1,'MinSeparation',3);
else
    local_tr = ref_raster_local<0 & islocalmin(ref_raster_local,1,'MinSeparation',3);
end
local_tr_ramp = polarity*vec([0:local_transient_window (local_transient_window-1):-1:0]); 
closest_local_tr = repmat(vec(local_tr_ramp),1,size(ref_raster_local,2),size(ref_raster_local,3));
% disp(size(closest_local_tr))
% disp(size(local_tr))
% disp(size(ref_raster))
closest_local_tr = local_tr.*(closest_local_tr+ref_raster_local);
closest_local_tr(closest_local_tr==0) = nan;

if polarity == 1
    [~,ref_tr_idx] = max(closest_local_tr,[],1);
    [~,closest_ref_tr_idx] = max(closest_local_tr,[],1,'linear');
else
    [~,ref_tr_idx] = min(closest_local_tr,[],1);
    [~,closest_ref_tr_idx] = min(closest_local_tr,[],1,'linear');
end
% get max values and idx   
ref_tr_mag = transpose(squeeze(ref_raster_local(closest_ref_tr_idx)));
ref_tr_idx = transpose(squeeze(ref_tr_idx));
% center it - this gives the position relative to the avg trg avg idx
rel_ref_tr_idx = ref_tr_idx-ceil(numel(local_tr_ramp)/2); 

% isolate the part of the corr_raster we want
corr_idx_i = repmat(vec(1:size(corr_raster,1)),1,size(corr_raster,2),size(corr_raster,3));
corr_idx_i = corr_idx_i(ismember(input_idx,corr_idx_of_int),:,:);
corr_idx_j = repmat(1:size(corr_raster,2),numel(corr_idx_of_int),1,size(corr_raster,3));
corr_idx_k = repmat(permute(1:size(corr_raster,3),[1 3 2]),numel(corr_idx_of_int),size(corr_raster,2),1);
if center_trial_transient == 0 % center the corr raster around the trg_avg_mean
    wind_idx_tp_centered = corr_idx_i + repmat(input_idx(avg_tr_idx),...
        numel(corr_idx_of_int),1,size(corr_raster,3));
    wind_idx_tp_centered(wind_idx_tp_centered<1 | wind_idx_tp_centered>numel(input_idx)) = nan;
else % center the corr raster around the trial-specific means if that's what we're doing
    wind_idx_tp_centered = corr_idx_i + repmat(permute(input_idx(rel_ref_tr_idx+avg_tr_idx),[3 2 1]),...
        numel(corr_idx_of_int),1,1);
    wind_idx_tp_centered(wind_idx_tp_centered<1 | wind_idx_tp_centered>numel(input_idx)) = nan;
end
corr_raster_centered = nan(size(wind_idx_tp_centered));
centered_ind = sub2ind(size(corr_raster),wind_idx_tp_centered(:),corr_idx_j(:),corr_idx_k(:));
corr_raster_centered(~isnan(centered_ind)) = corr_raster(centered_ind(~isnan(centered_ind)));

% now run the correlations
output.corr.corr_r = nan(size(corr_raster_centered,1),size(corr_raster_centered,2));
output.corr.p = nan(size(corr_raster_centered,1),size(corr_raster_centered,2));
for r = 1:size(ref_tr_mag,2)
    this_ref_tr_mag = ref_tr_mag(:,r);
    % we have to loop this way bc of the nan situation
    for t = 1:size(corr_raster_centered,1)
        this_corr_raster = squeeze(corr_raster_centered(t,r,:));
        keep_idx = ~isnan(this_ref_tr_mag) & ~isnan(this_corr_raster);
        if sum(keep_idx)>10
            [corr_r,p] = corr(this_ref_tr_mag(keep_idx),this_corr_raster(keep_idx,:));           
            output.corr.corr_r(t,r) = corr_r;
            output.corr.p(t,r) = p;
        end
    end
end


% results: dominant (largest magnitude) corr_r
[~,dominant_corr_r_idx] = max(abs(output.corr.corr_r),[],1);
[~,dominant_corr_r_lin_idx] = max(abs(output.corr.corr_r),[],1,'linear');
output.corr.dominant_corr_r = vec(output.corr.corr_r(dominant_corr_r_lin_idx));
output.corr.dominant_p = vec(output.corr.p(dominant_corr_r_lin_idx));
output.corr.dominant_corr_r_idx = vec(dominant_corr_r_idx);
output.corr.dominant_corr_r_input_idx = vec(corr_idx_of_int(dominant_corr_r_idx));


% results: max corr_r
[~,dominant_corr_r_idx] = max(output.corr.corr_r,[],1);
[~,dominant_corr_r_lin_idx] = max(output.corr.corr_r,[],1,'linear');
output.corr.max_corr_r = vec(output.corr.corr_r(dominant_corr_r_lin_idx));
output.corr.max_p = vec(output.corr.p(dominant_corr_r_lin_idx));
output.corr.max_corr_r_idx = vec(dominant_corr_r_idx);
output.corr.max_corr_r_input_idx = vec(corr_idx_of_int(dominant_corr_r_idx));


% results: min corr_r
[~,dominant_corr_r_idx] = min(output.corr.corr_r,[],1);
[~,dominant_corr_r_lin_idx] = min(output.corr.corr_r,[],1,'linear');
output.corr.min_corr_r = vec(output.corr.corr_r(dominant_corr_r_lin_idx));
output.corr.min_p = vec(output.corr.p(dominant_corr_r_lin_idx));
output.corr.min_corr_r_idx = vec(dominant_corr_r_idx);
output.corr.min_corr_r_input_idx = vec(corr_idx_of_int(dominant_corr_r_idx));


% some useful info
output.corr_raster.rel_idx = corr_idx_of_int;
output.corr_raster.raster_centered = corr_raster_centered;
output.corr_raster.mu = nanmean(corr_raster_centered,3);
output.corr_raster.sem = nanstd(corr_raster_centered,[],3)./sqrt(size(corr_raster_centered,3));
output.ref_raster.idx = ref_idx_of_int;
output.ref_raster.avg_tr_idx = avg_tr_idx;
output.ref_raster.trial_tr_idx = rel_ref_tr_idx; % relative to avg_tr_idx
output.ref_raster.mu = nanmean(ref_raster,3);
output.ref_raster.sem = nanstd(ref_raster,[],3)./sqrt(size(ref_raster,3));
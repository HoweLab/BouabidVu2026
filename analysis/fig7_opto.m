%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organization
addpath(fullfile(pwd,'common_functions'))
data_dir = 'D:';
mice = {'ADS6','ADS12','ADS13','ADS16','ADS17','ADS19'};
fib = cohort_fib_table(data_dir,mice);
% load corr hotspot
save_dir1 = fullfile(data_dir,'results','1_cross_corr');
corr_hotspot = load(fullfile(save_dir1,'cross_corr_dominant_results.mat'),'sig_moran','str');
% directory for saving (interim) results
save_dir7 = fullfile(data_dir,'results','7_opto');
if ~exist(save_dir7,'dir')
    mkdir(save_dir7)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. cross-correlation replication
% note: lags = ACh lag w/respect to DA, i.e. (+)lag means DA leads
sr = 18;
lag = 18;
concat_data = get_concat_data(mice,fib,data_dir,5*sr);
cc_results = get_cross_corr_results(concat_data,mice,fib,lag);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. cross-correlation replication maps
% size & striatum mask
voxel_size = 0.05;
str = get_striatum_vol_mask(...
    [min(fib.fiber_bottom_AP) max(fib.fiber_bottom_AP)],... % AP_range
    [min(fib.fiber_bottom_ML) max(fib.fiber_bottom_ML)],... % ML_range
    [min(fib.fiber_bottom_DV) max(fib.fiber_bottom_DV)],... % DV_range
    voxel_size);

% smooth maps and moran calculations
results = struct;

% vals
results.vals.r = cc_results.r.dominant;
results.vals.lag = cc_results.lag.dominant;
results.vals.cluster = vec(kmeans(...
    [results.vals.r,results.vals.lag],2));

% smooth maps
tmp = smooth_activity_map_interp_smoothed(results.vals.r, fib,...
    str.info.voxel_size,'AP_range',[min(str.info.AP) max(str.info.AP)],...
    'ML_range',[min(str.info.ML) max(str.info.ML)],...
    'DV_range',[min(str.info.DV) max(str.info.DV)],'gaussian_sigma',2);  
results.smooth = tmp.vol_01.smooth;

% moran
results.moran = local_morans_I(results.smooth,'weight_matrix',ones(21,21,21));

% null moran (this can take a while -- easier to run this on the side
% and save out results and then come back to it)
null_moran = local_morans_I_bootstrap_null(...
    fib,results.vals.r,voxel_size,...
    'AP_range',[min(str.info.AP) max(str.info.AP)],...
    'ML_range',[min(str.info.ML) max(str.info.ML)],...
    'DV_range',[min(str.info.DV) max(str.info.DV)]);


% significant hotspot
results.sig_moran = local_morans_I_sig_moran(results.smooth,...
    results.moran,null_moran,'dominant','mask',str.striatum_mask);

% correlation hotspot comparison
results.sig_moran.corr_hotspot_overlap = hotspot_comparison(...
    results.sig_moran.rand,results.sig_moran.vox,...
    str.info.DV,str.info.ML,str.info.AP,...
    corr_hotspot.sig_moran.rand,corr_hotspot.sig_moran.vox,...
    corr_hotspot.str.info.DV,corr_hotspot.str.info.ML,corr_hotspot.str.info.AP);

% save for convenience
results.str = str;
save(fullfile(save_dir7,'cross_corr_dominant_results.mat'),'-struct','results')
clear results

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. stim effects
stim_ta = get_stim_trg_avg(concat_data,sr); % useful to save out

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. stim effect maps
eta_idx = -9:45;
idx_of_int = -2:36;
opto_corr_hotspot = load(fullfile(save_dir7,'cross_corr_dominant_results.mat'),'sig_moran','str');
neuromods = {'ACh','DA'};
for n = 1:numel(neuromods)
    results = struct;
    
    % vals
    all_mags = zeros(size(fib,1),2);
    all_sig = zeros(size(fib,1),2);
    for m = 1:numel(mice)
        mouse = mice{m};
        mouse_raster = stim_ta.(mouse).(neuromods{n})(ismember(eta_idx,idx_of_int),:,:);
        max_mag = mean(max(mouse_raster),3);
        min_mag = mean(min(mouse_raster),3);
        all_mags(ismember(fib.mouse,mouse),:) = [vec(max_mag) vec(min_mag)];
        all_sig(ismember(fib.mouse,mouse),:) = abs(stim_ta.(mouse).([neuromods{n} '_sig']));
    end
    
    if strcmp(neuromods{n},'ACh') % we care about dip
        results.vals.mag = all_mags(:,2);
        results.vals.sig = all_sig(:,2);
    else % DA we care about peak
        results.vals.mag = all_mags(:,1);
        results.vals.sig = all_sig(:,1);
    end
    
    % smooth maps
    tmp = smooth_activity_map_interp_smoothed(results.vals.mag, fib,...
        str.info.voxel_size,'AP_range',[min(str.info.AP) max(str.info.AP)],...
        'ML_range',[min(str.info.ML) max(str.info.ML)],...
        'DV_range',[min(str.info.DV) max(str.info.DV)],'gaussian_sigma',2);  
    results.smooth = tmp.vol_01.smooth;

    % moran
    results.moran = local_morans_I(results.smooth,'weight_matrix',ones(21,21,21));
    
    % null moran (this can take a while -- easier to run this on the side
    % and save out results and then come back to it)
    null_moran = local_morans_I_bootstrap_null(...
        fib,results.vals.r,voxel_size,...
        'AP_range',[min(str.info.AP) max(str.info.AP)],...
        'ML_range',[min(str.info.ML) max(str.info.ML)],...
        'DV_range',[min(str.info.DV) max(str.info.DV)]);

    % significant hotspot
    results.sig_moran = local_morans_I_sig_moran(results.smooth,...
        results.moran,null_moran,'dominant','mask',str.striatum_mask);  
    
    % correlation hotspot comparison
    results.sig_moran.cohort_corr_hotspot_overlap = hotspot_comparison(...
        results.sig_moran.rand,results.sig_moran.vox,...
        str.info.DV,str.info.ML,str.info.AP,...
        opto_corr_hotspot.sig_moran.rand,opto_corr_hotspot.sig_moran.vox,...
        opto_corr_hotspot.str.info.DV,opto_corr_hotspot.str.info.ML,opto_corr_hotspot.str.info.AP);
    
    % save
    results.str = str;
    save(fullfile(save_dir7,[neuromods{n} '_stim_effect_results.mat']),'-struct','results')
end

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. results figs: ACh

results = load(fullfile(save_dir7,'ACh_stim_effect_results.mat'));
opto_corr_hotspot = load(fullfile(save_dir7,'cross_corr_dominant_results.mat'),'sig_moran','str');

% smooth maps
stim_outlines = get_mask_projection_outlines(results.sig_moran.map,...
    results.str,'apply_str_mask',1,'proj_orientations',{'axial','sagittal'});
corr_outlines = get_mask_projection_outlines(opto_corr_hotspot.sig_moran.map,...
    opto_corr_hotspot.str,'apply_str_mask',1,'proj_orientations',{'axial','sagittal'});
outlines.axial = [stim_outlines.axial; corr_outlines.axial];
outlines.sagittal = [stim_outlines.sagittal; corr_outlines.sagittal];
plot_smooth_maps(results.smooth,results.str,'outlines',outlines);

% in-v-out violin
plot_violin_in_out(results.vals.mag,fib,opto_corr_hotspot.str,opto_corr_hotspot.sig_moran.map)

% venn diagrams of hotspot comparisons (title has #voxels)
plot_hotspot_comparison_venn(results.sig_moran.cohort_corr_hotspot_overlap)

% sig pie
plot_sig_pie(results.vals.sig,[0 0 1])



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot_sig_pie
function plot_sig_pie(sig,sig_color)

% get counts
sig_count = sum(sig==1);
insig_count = sum(sig==0);
pie_vals = [sig_count insig_count];
pie_colors = [sig_color;.5 .5 .5];

figure
p = pie(pie_vals);
for i = 1:numel(pie_vals)
    p((i-1)*2+1).FaceColor = pie_colors(i,:);
    if pie_vals(i) >0
        perc_str = p(2*i).String;
        p(i*2).String = [perc_str ' (' num2str(pie_vals(i)) ')'];
    else
        p(i*2).String = '';
    end        
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_concat_data
function concat_data = get_concat_data(mice,fib,data_dir,nanpad)
concat_data = struct;
behav_fields = {'timestamp','reward','lick','stimDriver'};

for m = 1:numel(mice)
    mouse = mice{m};           
    str_rois = fib.ROI_orig(strcmp(fib.mouse,mouse));

    % preallocate
    concat_data.(mouse).date = [];
    concat_data.(mouse).DA = [];
    concat_data.(mouse).ACh = [];
    concat_data.(mouse).time_since_rew_del = [];
    concat_data.(mouse).time_since_rew_con = [];
    concat_data.(mouse).time_since_last_stim_on = [];
    for f = 1:numel(behav_fields)
        concat_data.(mouse).behav_DA.(behav_fields{f}) = [];
        concat_data.(mouse).behav_ACh.(behav_fields{f}) = [];
    end


    % loop over data, keeping only spontaneous movement sessions, i.e.,
    % non-reward sessions (non-reward session folders end in 'r')
    exp_dirs = dir(fullfile(data_dir,mouse));
    is_dirs = [exp_dirs.isdir];
    exp_dirs = {exp_dirs.name}';
    exp_dirs = exp_dirs(is_dirs);
    exp_dirs = exp_dirs(~startsWith(exp_dirs,'.'));
    for d = 1:numel(exp_dirs)
        data = load(fullfile(data_dir,mouse,exp_dirs{d},[mouse '_' exp_dirs{d} '.mat']));
        if ~isempty(data.ACh.Fc_exp) && ~isempty(data.DA.Fc_exp)
            % append
            concat_data.(mouse).DA = [concat_data.(mouse).DA; ...
                nan(nanpad,numel(str_rois)); data.DA.Fc_exp(:,str_rois)];
            concat_data.(mouse).ACh = [concat_data.(mouse).ACh;...
                nan(nanpad,numel(str_rois)); data.ACh.Fc_exp(:,str_rois)];
            concat_data.(mouse).date = [concat_data.(mouse).date;...
                repmat({'n/a'},nanpad,1); ...
                repmat(exp_dirs(d),size(data.DA.Fc_exp,1),1)];
            % time related to reward and stim
            rew_info = get_reward_info(data.behav_ACh);
            concat_data.(mouse).time_since_rew_del =...
                [concat_data.(mouse).time_since_rew_del;...
                nan(nanpad,1);rew_info.time_since_rew_del];
            concat_data.(mouse).time_since_rew_con =...
                [concat_data.(mouse).time_since_rew_con;...
                nan(nanpad,1);rew_info.time_since_rew_con];
            % stim
            stim_on = strfind(transpose(data.behav_ACh.stimDriver>0.1),[0 1]);
            time_since_last_stim_on = nan(size(data.behav_ACh.stimDriver));
            if ~isempty(stim_on)
                for st = 1:numel(stim_on)
                    s1 = stim_on(st);                                        
                    if st < numel(stim_on)
                        s2 = stim_on(st+1);
                    else
                        s2 = numel(time_since_last_stim_on);
                    end
                    time_since_last_stim_on(s1:(s2-1)) = ...
                        data.behav_ACh.timestamp(s1:(s2-1))-data.behav_ACh.timestamp(s1);
                end
            end
            concat_data.(mouse).time_since_last_stim_on =...
                [concat_data.(mouse).time_since_last_stim_on;...
                nan(nanpad,1);time_since_last_stim_on];
            % behav: useful for later
            for f = 1:numel(behav_fields)
                concat_data.(mouse).behav_DA.(behav_fields{f}) = [...
                    concat_data.(mouse).behav_DA.(behav_fields{f});...
                    nan(nanpad,1);data.behav_DA.(behav_fields{f})];
                concat_data.(mouse).behav_ACh.(behav_fields{f}) = [...
                    concat_data.(mouse).behav_ACh.(behav_fields{f});...
                    nan(nanpad,1);data.behav_ACh.(behav_fields{f})];
            end           
        end   
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_cross_corr_results
function cc_results = get_cross_corr_results(concat_data,mice,fib,lag,sr,varargin)
%%%  parse optional inputs %%%
ip = inputParser;
ip.addParameter('rew_ignore',6)
ip.addParameter('stim_ignore',3)
ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end

cc_results = struct;
cc_results.corr.r = nan(size(fib,1),numel(-lag:lag));
cc_results.corr.p = nan(size(fib,1),numel(-lag:lag));
cc_results.corr.n = nan(size(fib,1),numel(-lag:lag));

which_r = {'min','max','dominant'};
for i = 1:numel(which_r)
    cc_results.r.(which_r{i}) = nan(size(fib,1),1);
    cc_results.lag.(which_r{i}) = nan(size(fib,1),1);
    cc_results.p.(which_r{i}) = nan(size(fib,1),1);
    cc_results.n.(which_r{i}) = nan(size(fib,1),1);
end
for m = 1:numel(mice)
    mouse = mice{m};
    mouse_idx = find(ismember(fib.mouse,mouse));
    for r = 1:size(concat_data.(mouse).DA,2)
        this_DA = concat_data.(mouse).DA(:,r);
        this_ACh = concat_data.(mouse).ACh(:,r);
        % ignore stimulation and reward timepoints
        nan_filt = concat_data.(mouse).time_since_rew_del < 6 & ...
            concat_data.(mouse).time_since_rew_con < 6 & ...
            concat_data.(mouse).time_since_last_stim_on < 3;        
        this_DA(nan_filt) = nan;        
        this_ACh(nan_filt) = nan;
        for j = -lag:lag
            lag_ACh = [nan(-min(j,0),1);...
                this_ACh((max(j,0)+1):(end+min(j,0)),1);...
                nan(max(j,0),1)];
            keep_idx = ~isnan(this_DA) & ~isnan(lag_ACh);
            [corr_r,p] = corr(this_DA(keep_idx),lag_ACh(keep_idx));
            cc_results.corr.r(mouse_idx(r),j+lag+1) = corr_r;
            cc_results.corr.p(mouse_idx(r),j+lag+1) = p;
            cc_results.corr.n(mouse_idx(r),j+lag+1) = sum(keep_idx);
            disp([mouse ' ' num2str(r) ' ' num2str(j)])
        end
    end
end

% now get min, max, dominant r (and latencies)
% min
[min_r,min_lag] = min(cc_results.corr.r,[],2);
[~,min_lag_lin_idx] = min(cc_results.corr.r,[],2,'linear');
cc_results.r.min = min_r;
cc_results.lag.min = min_lag-lag-1;
cc_results.n.min = cc_results.corr.n(min_lag_lin_idx);
cc_results.p.min = cc_results.corr.p(min_lag_lin_idx);
% max
[max_r,max_lag] = max(cc_results.corr.r,[],2);
[~,max_lag_lin_idx] = max(cc_results.corr.r,[],2,'linear');
cc_results.r.max = max_r;
cc_results.lag.max = max_lag-lag-1;
cc_results.n.max = cc_results.corr.n(max_lag_lin_idx);
cc_results.p.max = cc_results.corr.p(max_lag_lin_idx);
% dominant 
cc_results.r.dominant = cc_results.r.max;
cc_results.lag.dominant = cc_results.lag.max;
cc_results.n.dominant = cc_results.n.max;
cc_results.p.dominant = cc_results.p.max;
min_better = abs(cc_results.r.min)>abs(cc_results.r.max);
cc_results.r.dominant(min_better) = cc_results.r.min(min_better);
cc_results.lag.dominant(min_better) = cc_results.lag.min(min_better);
cc_results.n.dominant(min_better) = cc_results.n.min(min_better);
cc_results.p.dominant(min_better) = cc_results.p.min(min_better);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_stim_trg_avg
function output = get_stim_trg_avg(concat_data,sr,varargin)

%%%  parse optional inputs %%%
ip = inputParser;
ip.addParameter('eta_idx',-9:45)
ip.addParameter('idx_of_int',-2:36)
ip.addParameter('n_rew_excl_s',5)
ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end

neuromods = {'DA','ACh'};
mice = fieldnames(concat_data);
for m = 1:numel(mice)
    mouse = mice{m};
    for n = 1:numel(neuromods)
        rew_info = get_reward_info(concat_data.(mouse).(['behav_' neuromods{n}]));
        rew_stim = eventTriggeredAverage(double(concat_data.(mouse).(['behav_' neuromods{n}]).stimDriver>0.1),rew_info.rew_del,1,18,'nullDistr',0);
        rew_stim = sum(squeeze(rew_stim.activity))>0;
        excl{n} = [vec(repmat(0:(n_rew_excl_s*sr),numel(rew_info.rew_del),1)+rew_info.rew_del);... % excl 5s of reward-related activity
            vec(repmat(0:(n_rew_excl_s*sr),sum(~isnan(rew_info.rew_lick)),1)+rew_info.rew_lick(~isnan(rew_info.rew_lick)))]; % excl 5s of reward-related activity
        stim_on = strfind(transpose(double(concat_data.(mouse).(['behav_' neuromods{n}]).stimDriver>0.1)),[zeros(1,100) 1])+100;
        excl{n} = unique([excl{n}; vec(repmat(1:18,numel(stim_on),1) + vec(stim_on))]);
        iti_stim = eventTriggeredAverage(concat_data.(mouse).(['behav_' neuromods{n}]).reward,stim_on,-18,18,'nullDistr',0);
        iti_stim = sum(squeeze(iti_stim.activity))==0;
        stim_events{n} = stim_on(iti_stim==1);
    end

    % now check numbers
    while numel(stim_events{1}) ~= numel(stim_events{2})
        % which one is missing a channel
        [~,ch_m] = min(cellfun(@(x) numel(x),stim_events));
        % let's figure out what the "offset" is between
        % this and the other channel
        [~,closest_idx] = arrayfun(@(x) min(abs(x-stim_events{2/ch_m})),stim_events{ch_m});
        idx_offset = mode(stim_events{ch_m}-stim_events{2/ch_m}(closest_idx));
        % now find the missing one, and add that in there
        % with the offset
        missing_idx = setdiff(stim_events{2/ch_m},...
            [stim_events{ch_m}; ...
            stim_events{ch_m}+1; ...
            stim_events{ch_m}-1]);
        stim_events{ch_m} = sort([stim_events{ch_m}; missing_idx+idx_offset]);
    end
    
    % event-triggered average
    for n = 1:numel(neuromods)                    
        eta = eventTriggeredAverage(concat_data.(mouse).(neuromods{n}),...
            stim_events{n},eta_idx(1),eta_idx(end),'nullDistr',1,'bootstrapN',10000,'otherExcl',excl{n});
        eta = eta_sig(eta,'sig_idx',[find(eta_idx==idx_of_int(1)) find(eta_idx==idx_of_int(end))]);
        output.(mouse).(neuromods{n}) = eta.activity;  
        output.(mouse).([neuromods{n} '_sig']) = eta.sig_cols;  
    end
end
end





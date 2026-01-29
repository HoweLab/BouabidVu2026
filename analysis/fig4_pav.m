%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organization
addpath(fullfile(pwd,'common_functions'))
data_dir = 'D:';
mice = {'AD1','AD2','AD3'};
fib = cohort_fib_table(data_dir,mice);
save_dir1 = fullfile(data_dir,'results','1_cross_corr');
corr_hotspot = load(fullfile(save_dir1,'cross_corr_dominant_results.mat'),'sig_moran','str');
save_dir4 = fullfile(data_dir,'results','4_pav');
if ~exist(save_dir4,'dir')
    mkdir(save_dir4)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. get concatenated pav data 
epochs = {'pav','led';...
    'led','led';...
    'pav','tone';...
    'tone','tone'};
pav_data = get_all_pav_cue_data(mice,fib,data_dir,epochs); % might be useful to save this


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. cross correlation of trial-by-trial ACh w trial-by-trial DA transient
% interested in peaks for pav, dips for extinction;
DA_tr = double(strcmp(epochs(:,1),'pav')) + -1*double(~strcmp(epochs(:,1),'pav'));
pav_cc = get_all_pav_cc(pav_data,fib,epochs,DA_tr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. maps

% size & striatum mask
voxel_size = 0.05;
str = get_striatum_vol_mask(...
    [min(fib.fiber_bottom_AP) max(fib.fiber_bottom_AP)],... % AP_range
    [min(fib.fiber_bottom_ML) max(fib.fiber_bottom_ML)],... % ML_range
    [min(fib.fiber_bottom_DV) max(fib.fiber_bottom_DV)],... % DV_range
    voxel_size);

for e = 1:size(epochs,1)

    % smooth maps and moran calculations
    results = struct;

    % vals
    results.vals.r = pav_cc.dominant.(strjoin(epochs(e,1:2),'_')).corr_r;
    results.vals.p = pav_cc.dominant.(strjoin(epochs(e,1:2),'_')).p;
    results.vals.lat = pav_cc.dominant.(strjoin(epochs(e,1:2),'_')).lat;
    results.vals.dff_ach =  pav_cc.dominant.(strjoin(epochs(e,1:2),'_')).dff_ach;
    results.vals.sig_da =  pav_cc.dominant.(strjoin(epochs(e,1:2),'_')).sig_da;
    results.vals.cluster = nan(size(results.vals.sig_da));
    results.vals.cluster(results.vals.sig_da==1) = ...
        vec(kmeans([results.vals.r(results.vals.sig_da==1),results.vals.lat(results.vals.sig_da==1)],2));

    % smooth maps
    these_vals = results.vals.r(results.vals.sig_da==1);
    this_fib = fib(results.vals.sig_da==1,:);
    this_fib.ROI = vec(1:size(this_fib,1));
    tmp = smooth_activity_map_interp_smoothed(these_vals, this_fib,...
        str.info.voxel_size,'AP_range',[min(str.info.AP) max(str.info.AP)],...
        'ML_range',[min(str.info.ML) max(str.info.ML)],...
        'DV_range',[min(str.info.DV) max(str.info.DV)],'gaussian_sigma',2);  
    results.smooth = tmp.vol_01.smooth;

    % moran
    results.moran = local_morans_I(results.smooth,'weight_matrix',ones(21,21,21));

    % null moran (this can take a while -- easier to run this on the side
    % and save out results and then come back to it)
    null_moran = local_morans_I_bootstrap_null(...
        this_fib,these_vals,voxel_size,...
        'AP_range',[min(str.info.AP) max(str.info.AP)],...
        'ML_range',[min(str.info.ML) max(str.info.ML)],...
        'DV_range',[min(str.info.DV) max(str.info.DV)]);


    % significant hotspot
    results.sig_moran = local_morans_I_sig_moran(results.smooth,...
        results.moran,null_moran,'dominant','mask',str.striatum_mask);  

    % save for convenience
    results.str = str;
    save(fullfile(save_dir4,[strjoin(epochs(e,1:2),'_') '_dominant_results.mat']),'-struct','results')
end

%%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. results figs
for e = 1:2 % plot light trials results for paper
    results = load(fullfile(save_dir4,[strjoin(epochs(e,1:2),'_') '_dominant_results.mat']));
    this_fib = fib(results.vals.sig_da==1,:);
    this_fib.ROI = vec(1:size(this_fib,1));
    
    % scatter plot
    plot_lag_corr_scatter(results.vals.r(results.vals.sig_da==1),...
        results.vals.lat(results.vals.sig_da==1)/18*1000,...
        'lat_bins',[-500:(1000/18):500],'clust_id',...
        results.vals.cluster(results.vals.sig_da==1));

    % smooth maps
    ur_outlines = get_mask_projection_outlines(results.sig_moran.map,...
        results.str,'apply_str_mask',1,'proj_orientations',{'axial','sagittal'});
    corr_outlines = get_mask_projection_outlines(corr_hotspot.sig_moran.map,...
        corr_hotspot.str,'apply_str_mask',1,'proj_orientations',{'axial','sagittal'});
    outlines.axial = [ur_outlines.axial; corr_outlines.axial];
    outlines.sagittal = [ur_outlines.sagittal; corr_outlines.sagittal];
    plot_smooth_maps(results.smooth,results.str,'outlines',outlines);

    % in-v-out violin
    plot_violin_in_out(results.vals.r(results.vals.sig_da==1),...
        this_fib,corr_hotspot.str,corr_hotspot.sig_moran.map)
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_all_pav_cue_data
% epochs = nx2 cell array where first column is task (e.g., 'pav'), 
% 2nd column is trial type, e.g., 'led'
% notes: 
% 1. pav = learning phase; tone = tone extinction; led = led extinction
function pav_data = get_all_pav_cue_data(mice,fib,data_dir,epochs,varargin)

%%%  parse optional inputs %%%
ip = inputParser;
ip.addParameter('n_days',3); % consider last n_days
ip.addParameter('eta_idx',-9:54)
ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end

% hard-coded info
cues = {'led','tone'}; % cue1 = led, cue2 = tone
cue_offset = [0 7]; % there's a 7-frame mechanical tone onset delay
neuromods = {'ACh','DA'};


pav_data = struct;
for m = 1:numel(mice)
    mouse = mice{m};
    str_rois = fib.ROI_orig(strcmp(fib.mouse,mouse));

    % loop over data
    exp_dirs = dir(fullfile(data_dir,mouse));
    is_dirs = [exp_dirs.isdir];
    exp_dirs = {exp_dirs.name}';
    exp_dirs = exp_dirs(is_dirs);
    exp_dirs = exp_dirs(~startsWith(exp_dirs,'.'));
    
    for e = 1:size(epochs,1)
        cue_c = find(strcmp(cues,epochs{e,2}));
        
        % get relevant recordings
        task_dirs = [];
        for d = 1:numel(exp_dirs)
            % check task first
            data = load(fullfile(data_dir,mouse,exp_dirs{d},[mouse '_' exp_dirs{d} '.mat']),'task');
            if strcmp(data.task,epochs{e,1}) % if the right task, add to list
                task_dirs = [task_dirs;exp_dirs(d)];
            end
        end
        task_dirs = task_dirs((end-n_days+1):end);
        
        % tmp cells for triggered averages
        tmp = struct;
        tmp.DA = cell(0,0);
        tmp.ACh = cell(0,0);
        % loop over relevant recordings
        for d = 1:numel(task_dirs)
            data = load(fullfile(data_dir,mouse,exp_dirs{d},[mouse '_' exp_dirs{d} '.mat']));
            if ~isempty(data.DA.Fc_exp_hp) && ~isempty(data.ACh.Fc_exp_hp)
                for n = 1:numel(neuromods)
                    cue_info = pav2cue_getTrialTimes(data.(['behav_' neuromods{n}]));
                    cue_events = cue_info.cueStart(cue_info.cue==cue_c) + cue_offset(cue_c);                    
                    other_excl = pav2cue_unrewarded_ITI(data.(['behav_' neuromods{n}]));
                    other_excl = find(other_excl==0);

                    eta = eventTriggeredAverage(data.(neuromods{n}).Fc_exp_hp(:,str_rois),...
                        cue_events,eta_idx(1),eta_idx(end),'nullDistr',1,'bootstrapN',1,...
                        'otherExcl',other_excl);
                    tmp.(neuromods{n})(end+1,:) = [{eta} {data.(neuromods{n}).Fc_exp_hp(:,str_rois)}];
                end
            end
        end
        for n = 1:numel(neuromods)
            eta = eta_combine(tmp.(neuromods{n}));
            eta = eta_sig(eta,'sig_idx',[find(eta_idx==0) find(eta_idx==18)]); % significance assessed 0-1s
            pav_data.(mouse).(strjoin(epochs(e,1:2),'_')).(neuromods{n}) = eta;
        end
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_all_pav_cc
function pav_cc = get_all_pav_cc(pav_data,fib,epochs,DA_tr,varargin)

%%%  parse optional inputs %%%
ip = inputParser;
ip.addParameter('eta_idx',-9:54)
ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end
mice = fieldnames(pav_data);
which_r = {'min','max','dominant'};

% initialize
for i = 1:numel(which_r)
    for e = 1:size(epochs,1)
        pav_cc.(which_r{i}).(strjoin(epochs(e,1:2),'_')).corr_r = nan(size(fib,1),1);
        pav_cc.(which_r{i}).(strjoin(epochs(e,1:2),'_')).p = nan(size(fib,1),1);
        pav_cc.(which_r{i}).(strjoin(epochs(e,1:2),'_')).lat = nan(size(fib,1),1);
        pav_cc.(which_r{i}).(strjoin(epochs(e,1:2),'_')).dff_ach = nan(size(fib,1),1);
        pav_cc.(which_r{i}).(strjoin(epochs(e,1:2),'_')).sig_da = nan(size(fib,1),1);
    end
end

% loop
for m = 1:numel(mice)
    mouse = mice{m};   
    mouse_idx = ismember(fib.mouse,mouse);
    for e = 1:size(epochs,1)
        da_act = pav_data.(mouse).(strjoin(epochs(e,1:2),'_')).DA.activity;
        ach_act = pav_data.(mouse).(strjoin(epochs(e,1:2),'_')).ACh.activity;
        % if DA has the significant feature we're interested in
        if DA_tr(e) == 1 % if we're interested in DA peak
            idx_of_int = 0:18; % significance assessed 0-1s
            pav_data.(mouse).(strjoin(epochs(e,1:2),'_')).DA = ...
                eta_sig(pav_data.(mouse).(strjoin(epochs(e,1:2),'_')).DA,...
                'sig_idx',[find(eta_idx==idx_of_int(1)) find(eta_idx==idx_of_int(end))]);
            sig_da = abs(double(pav_data.(mouse).(strjoin(epochs(e,:),'_')).DA.sig_cols(:,1)));          

        else % significant dip
            idx_of_int = 0:27; % significance assessed 0-1.5s
            pav_data.(mouse).(strjoin(epochs(e,1:2),'_')).DA = ...
                eta_sig(pav_data.(mouse).(strjoin(epochs(e,1:2),'_')).DA,...
                'sig_idx',[find(eta_idx==idx_of_int(1)) find(eta_idx==idx_of_int(end))]); 
            sig_da = abs(double(pav_data.(mouse).(strjoin(epochs(e,:),'_')).DA.sig_cols(:,2)));
        end
        
        % correlate ACh with DA peak magnitude
        output = raster_transient_correlation(da_act,ach_act,1,...
            'input_idx',eta_idx,'ref_idx_of_int',idx_of_int,'corr_idx_of_int',-9:9,...
            'local_transient_window',9,'center_trial_transient',1); 
        
        % record
        for i = 1:numel(which_r)
            pav_cc.(which_r{i}).(strjoin(epochs(e,1:2),'_')).corr_r(mouse_idx) = ...
                output.corr.([which_r{i} '_corr_r']);
            pav_cc.(which_r{i}).(strjoin(epochs(e,1:2),'_')).p(mouse_idx) = ...
                output.corr.([which_r{i} '_p']);
            pav_cc.(which_r{i}).(strjoin(epochs(e,1:2),'_')).lat(mouse_idx) = ...
                output.corr.([which_r{i} '_corr_r_input_idx']);       
            pav_cc.(which_r{i}).(strjoin(epochs(e,1:2),'_')).dff_ach(mouse_idx) = ...
                output.corr_raster.mu(sub2ind(size(output.corr_raster.mu),...
                output.corr.([which_r{i} '_corr_r_idx']),...
                vec(1:size(output.corr_raster.mu,2))));
            pav_cc.(which_r{i}).(strjoin(epochs(e,1:2),'_')).sig_da(mouse_idx) = sig_da;
        end
    end
end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pav2cue_getTrialTimes
%
% function output = pav2cue_getTrialTimes(behav)
%
% This function takes as input the path to a behavior file, or the struct
% itself, and returns a struct, with fields:
%   -trial          the trial number
%   -cue            1 or 2 (see behav.experimentSetup.exp.cue_freq to know
%                   which is which)
%   -cue_ID         the frequency in kHz (-1 if LED)
%   -level_analog   brightness or relative volume (actual level)
%   -level_cat      brightness or relative volume (categorical, 
%                   coded as 1, 2, 3, least to most salient)
%   -cueStart       index for cue starts 
%   -cueEnd         index for cue ends
%   -rew            whether the trial was rewarded (1) or not (0)
%   -rewOn          index for reward delivery
%   -rewLick        index of first lick after reward delivery
%   -trialTimes     trial time in s (how long from start of cue to end of cue)
%   -unpredRewOn    onset index of unpredicted reward
%   -unpredRewLIck  index for first lick after unpredicted reward
%   -unpredRewSize  size of unpredicted rewards
%   -fr             framerate (NOTE: ASSUMES INTEGER RATE)
% 
%
% Mai-Anh, updated 12/10/2021
% updated 2/1/2022 to handle the case where the imaging camera starts after
%                      the task has already begun
% updated 2/22/22 for more general use
%
function output = pav2cue_getTrialTimes(behav)

% load if it's a path
if ~isstruct(behav)
    if ~endsWith(behav,'.mat')
        behav = [behav '.mat'];
    end
    behav = load(behav);
end
trialInfo = behav.experimentSetup.exp.trials;
rewInfo = getRewInfo(behav);

% preallocate
output.trial = trialInfo.trialNum;
output.cue = trialInfo.cueRL;
output.cue_ID = transpose(behav.experimentSetup.exp.cue_freq(trialInfo.cueRL));
output.level_analog = transpose(behav.experimentSetup.exp.cue_relativeVol(trialInfo.cueRL));
output.level_cat = nan(size(output.level_analog));
output.cueStart = nan(size(output.trial));
output.rew = trialInfo.rew;
output.rewOn = nan(size(output.trial));
output.rewLick = nan(size(output.trial));
output.cueEnd = nan(size(output.trial));

cueFields = {'',''};
for i = 1:2
    if behav.experimentSetup.exp.cue_freq(i)==-1
        if isfield(behav.experimentSetup.exp,'stim_driver_freq') &&...
                behav.experimentSetup.exp.stim_driver_freq >0
            cueFields{i} = 'stimulus_led';
        else
            cueFields{i} = 'stimulus_led_analog';
        end
    else
        cueFields{i} = ['stimulus_sound' num2str(i)];
    end
end

% cue onsets and offsets & IDs so that we can match in case the first trial
% isn't recorded
stim_on = [];
stim_off = [];
stim_on_id = [];
stim_off_id = [];
stim_thresh = .1;
for i = 1:2
    stim_field = behav.(cueFields{i});
    stim_on = [stim_on; find(diff(stim_field>stim_thresh)==1)+1];
    stim_on_id = [stim_on_id; ones(numel(find(diff(stim_field>stim_thresh)==1)),1)*i];
    stim_off = [stim_off; find(diff(stim_field>stim_thresh)==-1)+1];
    stim_off_id = [stim_off_id; ones(numel(find(diff(stim_field>stim_thresh)==-1)),1)*i];
end
% in case we have mismatch
[stim_on,on_sort_idx]= sort(stim_on);
stim_on_id = stim_on_id(on_sort_idx);
[stim_off,off_sort_idx] = sort(stim_off);
stim_off_id = stim_off_id(off_sort_idx);
if stim_off(1)<stim_on(1)
    stim_off = stim_off(2:end);
    stim_off_id = stim_off_id(2:end);
end
if stim_on(end)>stim_off(end)
    stim_on = stim_on(1:end-1);
    stim_on_id = stim_on_id(1:end-1);
end
startTrial = strfind(trialInfo.cueRL', stim_on_id');
output.cueStart(startTrial:(startTrial+numel(stim_on)-1),1) = stim_on;
output.cueEnd(startTrial:(startTrial+numel(stim_off)-1),1) = stim_off;
output.trialTimes(startTrial:(startTrial+numel(stim_off)-1),1) = behav.timestamp(stim_off)-behav.timestamp(stim_on);

% let's add the reward onsets
for i = 1:numel(output.trial)
    if sum(rewInfo.rewOn>output.cueStart(i) & rewInfo.rewOn<output.cueEnd(i)) == 1
        output.rewOn(i) = rewInfo.rewOn((rewInfo.rewOn>output.cueStart(i) & rewInfo.rewOn<output.cueEnd(i)));
        output.rewLick(i) = rewInfo.lickOnset((rewInfo.rewOn>output.cueStart(i) & rewInfo.rewOn<output.cueEnd(i)));
    end
end

% level category
for i = 1:2
    thisCue = output.cue==i;
    levels = sort(unique(output.level_analog(thisCue)));
    for j = 1:numel(levels)
        thisCueLevel = output.cue==i & output.level_analog == levels(j);
        output.level_cat(thisCueLevel) = j;
    end
end    

% truncate uncompleted trials
output_fields = fieldnames(output);
keep_trials = ~isnan(output.cueEnd);
for f = 1:numel(output_fields)
    output.(output_fields{f}) = output.(output_fields{f})(keep_trials);
end

% unpredicted rewards
% include a +/-1 allowance to figure out which rewards were unpredicted
[UR,idx] = setdiff(rewInfo.rewOn,output.rewOn);
URsize = rewInfo.rewSize(idx);
output.unpredRewOn = UR;
output.unpredRewLick = rewInfo.lickOnset(idx);
output.unpredRewSize = URsize;

% add frame rate (assume integer)
output.fr = round(1/nanmean(diff(behav.timestamp)));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pav2cue_unrewarded_ITI
% Mai-Anh Vu, updated 8/16/2024
function eligible_timepoints = pav2cue_unrewarded_ITI(behav)

% load behav if necessary
if ischar(behav)
    if ~endsWith(behav,'.mat')
        behav = [behav '.mat'];
    end
    behav = load(behav);
end

% get ITI periods
eligible_timepoints = zeros(size(behav.timestamp));
% to be extra sure, do this a few ways
stimulus_fields = fieldnames(behav);
stimulus_fields = stimulus_fields(contains(stimulus_fields,'stim') & ~contains(stimulus_fields,'count'));
if isfield(behav.experimentSetup.exp,'cue_relativeVol')
    stim_thresh = min(behav.experimentSetup.exp.cue_relativeVol)/2;
else
    stim_thresh = 0.5;
end
if ~isempty(stimulus_fields)
    for s = 1:numel(stimulus_fields)
        eligible_timepoints = eligible_timepoints + ...
            double(behav.(stimulus_fields{s})>stim_thresh);
    end
    eligible_timepoints = double(eligible_timepoints==0);
end

% take out the ITIs with unpred reward
rew_info = getRewInfo(behav);
rew_outside_cue = rew_info.rewOn(eligible_timepoints(rew_info.rewOn)==1);
for r = 1:numel(rew_outside_cue)
    this_ur = rew_outside_cue(r);
    if find(eligible_timepoints==0,1,'first') > this_ur
        iti_ur_on = 1;
    else
        iti_ur_on = find(eligible_timepoints(1:this_ur)==0,1,'last')+1;
    end
    if find(eligible_timepoints==0,1,'last') < this_ur
        iti_ur_off = numel(eligible_timepoints);
    else
        iti_ur_off = find(eligible_timepoints(this_ur:end)==0,1,'first') + this_ur - 2;
    end
    eligible_timepoints(iti_ur_on:iti_ur_off) = 0;
end
end

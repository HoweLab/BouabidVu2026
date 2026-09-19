%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organization
addpath(fullfile(pwd,'common_functions'))
data_dir = 'G:';
coh_mice.cohort1 =  {'UG27','UG28','UG29','UG30','UG31'};
coh_mice.cohort2 = {'AD1','AD2','AD3'};
coh_mice.cohort5 = {'609','610','813','816','875',}; % ignoring 319
mice = struct2cell(structfun(@(x) x(:),coh_mice,'UniformOutput',false));
mice = vertcat(mice{:});
fib = cohort_fib_table(data_dir,mice);
% load corr hotspot
save_dir1 = fullfile(data_dir,'results','1_cross_corr');
%corr_hotspot = 
% directory for saving (interim) results
save_dir3 = fullfile(data_dir,'results','3_unpred_rew');
if ~exist(save_dir3,'dir')
    mkdir(save_dir3)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. get concatenated reward data 
rew_data = load_if_exist(fullfile(save_dir3,'rew_data.mat'));
if isempty(fieldnames(rew_data))
    rew_data = get_all_rew_data(mice,fib,data_dir,'Fc_field','Fc',...
        'Fc_artifact_mask','artifact_mask'); % useful to save this
    save(fullfile(save_dir3,'rew_data.mat'),'-struct','rew_data','-v7.3')
end

 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. cross correlation of trial-by-trial ACh w trial-by-trial DA peak 
rew_da = load_if_exist(fullfile(save_dir3,'rew_DA_peak_corr.mat'));
if isempty(fieldnames(rew_da))
    rew_da = get_DA_peak_rew_cc(rew_data,10000);
    save(fullfile(save_dir3,'rew_DA_peak_corr.mat'),'-struct','rew_da')
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. cross correlation of trial-by-trial DA w trial-by-trial ACh peak 
rew_ach = load_if_exist(fullfile(save_dir3,'rew_ACh_peak_corr.mat'));
if isempty(fieldnames(rew_ach))
    rew_ach = get_ACh_peak_rew_cc(rew_data,10000);
    save(fullfile(save_dir3,'rew_ACh_peak_corr.mat'),'-struct','rew_ach')
end
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 3. maps
% 
% % size & striatum mask
% voxel_size = 0.05;
% str = get_striatum_vol_mask(...
%     [min(fib.fiber_bottom_AP) max(fib.fiber_bottom_AP)],... % AP_range
%     [min(fib.fiber_bottom_ML) max(fib.fiber_bottom_ML)],... % ML_range
%     [min(fib.fiber_bottom_DV) max(fib.fiber_bottom_DV)],... % DV_range
%     voxel_size);
% 
% % smooth maps and moran calculations
% results = struct;
% 
% % vals
% results.vals.r = rew_cc.dominant.corr_r;
% results.vals.p = rew_cc.dominant.p;
% results.vals.lat = rew_cc.dominant.lat;
% results.vals.dff_ach = rew_cc.dominant.dff_ach;
% results.vals.cluster = vec(kmeans([results.vals.r,results.vals.lat],2));
% 
% % smooth maps
% tmp = smooth_activity_map_interp_smoothed(results.vals.r, fib,...
%     str.info.voxel_size,'AP_range',[min(str.info.AP) max(str.info.AP)],...
%     'ML_range',[min(str.info.ML) max(str.info.ML)],...
%     'DV_range',[min(str.info.DV) max(str.info.DV)],'gaussian_sigma',2);  
% results.smooth = tmp.vol_01.smooth;
% 
% % moran
% results.moran = local_morans_I(results.smooth,'weight_matrix',ones(21,21,21));
% 
% % null moran (this can take a while -- easier to run this on the side
% % and save out results and then come back to it)
% null_moran = local_morans_I_bootstrap_null(...
%     fib,results.vals.r,voxel_size,...
%     'AP_range',[min(str.info.AP) max(str.info.AP)],...
%     'ML_range',[min(str.info.ML) max(str.info.ML)],...
%     'DV_range',[min(str.info.DV) max(str.info.DV)]);
% 
% % significant hotspot
% results.sig_moran = local_morans_I_sig_moran(results.smooth,...
%     results.moran,null_moran,'dominant','mask',str.striatum_mask);
% 
% % correlation hotspot comparison
% results.sig_moran.corr_hotspot_overlap = hotspot_comparison(...
%     results.sig_moran.rand,results.sig_moran.vox,...
%     str.info.DV,str.info.ML,str.info.AP,...
%     corr_hotspot.sig_moran.rand,corr_hotspot.sig_moran.vox,...
%     corr_hotspot.str.info.DV,corr_hotspot.str.info.ML,corr_hotspot.str.info.AP);
% 
% % save for convenience
% results.str = str;
% save(fullfile(save_dir3,'unpred_rew_dominant_results.mat'),'-struct','results')
%     
% %%   
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 4. results figs
% results = load(fullfile(save_dir3,'unpred_rew_dominant_results.mat'));
% 
% % scatter plot
% plot_lag_corr_scatter(results.vals.r,results.vals.lat/18*1000,...
%     'lat_bins',[-500:(1000/18):500],'clust_id',results.vals.cluster);
% 
% % smooth maps
% ur_outlines = get_mask_projection_outlines(results.sig_moran.map,...
%     results.str,'apply_str_mask',1,'proj_orientations',{'axial','sagittal'});
% corr_outlines = get_mask_projection_outlines(corr_hotspot.sig_moran.map,...
%     corr_hotspot.str,'apply_str_mask',1,'proj_orientations',{'axial','sagittal'});
% outlines.axial = [ur_outlines.axial; corr_outlines.axial];
% outlines.sagittal = [ur_outlines.sagittal; corr_outlines.sagittal];
% plot_smooth_maps(results.smooth,results.str,'outlines',outlines);
% 
% % in-v-out violin
% plot_violin_in_out(results.vals.r,fib,corr_hotspot.str,corr_hotspot.sig_moran.map)
% 
% % venn diagrams of hotspot comparisons (title has #voxels)
% plot_hotspot_comparison_venn(results.sig_moran.corr_hotspot_overlap)
    

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_all_rew_data
function rew_data = get_all_rew_data(mice,fib,data_dir,varargin)

    %%%  parse optional inputs %%%
    ip = inputParser;
    ip.addParameter('eta_idx',-18:27);
    ip.addParameter('channel_names',{'ACh','DA'});  % channel names
    ip.addParameter('Fc_field','Fc');               % which DF/F field to use
    ip.addParameter('Fc_artifact_mask',[]);         % the field which contains the artifact mask; leave blank if N/A
    ip.addParameter('datafile_suffix','');          % datafile format [mouse]_[expdir][datafile_suffix].mat; default '', assume format MOUSE_EXPDIR.mat 
    ip.parse(varargin{:});
    for j=fields(ip.Results)'
        eval([j{1} '=ip.Results.' j{1} ';']);
    end

    rew_data = struct;
    for m = 1:numel(mice)
        mouse = mice{m};    
        mouse_field = get_mouse_field(mouse);
        str_rois = fib.ROI_orig(strcmp(fib.mouse,mouse));

        % loop over data, keeping only unpred rew sessions, i.e.,
        % folders end in 'r'
        exp_dirs = dir(fullfile(data_dir,mouse));
        is_dirs = [exp_dirs.isdir];
        exp_dirs = {exp_dirs.name}';
        exp_dirs = exp_dirs(is_dirs);
        exp_dirs = exp_dirs(~startsWith(exp_dirs,'.'));
        
        % further filter depending on the 
        if sum(endsWith(exp_dirs,'r')) > 0
            exp_dirs = exp_dirs(endsWith(exp_dirs,'r'));
            within_task = 0;
        else
            within_task = 1; % if there aren't 'r' experiments, then the unpredicted rewards are embedded within the Pav task
        end
        
        % tmp cells for reward triggered averages
        tmp = struct;
        for n = 1:2
            tmp.(channel_names{n}) = cell(0,0);
        end        
        for d = 1:numel(exp_dirs)
            data = load(fullfile(data_dir,mouse,exp_dirs{d},[mouse '_' exp_dirs{d} datafile_suffix '.mat']));
            if ~isempty(data.(channel_names{1}).(Fc_field)) && ~isempty(data.(channel_names{2}).(Fc_field))
                for n = 1:numel(channel_names)
                    if within_task == 0
                        rew_info = get_reward_info(data.(['behav_' channel_names{n}]));
                        these_rew = rew_info.rew_lick(~isnan(rew_info.rew_lick)); % trigger on consumption                        
                    else
                        rew_info = pav2cue_getTrialTimes(data.(['behav_' channel_names{n}]));
                        these_rew = rew_info.unpredRewLick(~isnan(rew_info.unpredRewLick)); % trigger on consumption
                    end
                    this_fc = data.(channel_names{n}).(Fc_field);
                    if ~isempty(Fc_artifact_mask)
                        this_fc(data.(channel_names{n}).(Fc_artifact_mask)) = nan;                        
                    end
                    this_fc = this_fc(:,str_rois);
                    eta = eventTriggeredAverage(this_fc, these_rew, eta_idx(1),eta_idx(end),'nullDistr',1,'bootstrapN',1);
                    tmp.(channel_names{n})(end+1,:) = [{eta} {data.(channel_names{n}).(Fc_field)(:,str_rois)}];
                end
            end   
        end
        for n = 1:numel(channel_names)
            eta = eta_combine(tmp.(channel_names{n}));
            eta = eta_sig(eta,'sig_idx',[find(eta_idx==0) find(eta_idx==18)]); % significance assessed 0-1s
            rew_data.(mouse_field).(channel_names{n}) = eta;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_DA_peak_rew_cc(
function rew_da = get_DA_peak_rew_cc(rew_data,null_it,varargin)

    %%%  parse optional inputs %%%
    ip = inputParser;
    ip.addParameter('eta_idx',-18:27); 
    ip.addParameter('corr_idx',0:9); 
    ip.addParameter('ref_idx_of_int',0:18); 
    ip.parse(varargin{:});
    for j=fields(ip.Results)'
        eval([j{1} '=ip.Results.' j{1} ';']);
    end
    
    % loop
    mice = fieldnames(rew_data);

    for m = 1:numel(mice)
        mouse = mice{m};         
        disp(mouse)
    %     mouse_idx = ismember(fib.mouse,mouse);
        da_act = rew_data.(mouse).DA.activity;
        ach_act = rew_data.(mouse).ACh.activity;
        % correlate ACh with DA peak magnitude
        mouse_corr = raster_transient_correlation(da_act,ach_act,1,...
            'input_idx',eta_idx,'ref_idx_of_int',ref_idx_of_int,...
                'corr_idx_of_int',corr_idx,'local_transient_window',9);      
        null_r = nan(null_it,size(da_act,2));
        % now null: shuffle the trials
        for i = 1:null_it      
            if rem(i,500) == 0
                disp(['     ' num2str(i)])
            end
            this_ach = ach_act(:,:,randperm(size(ach_act,3)));
            this_corr = raster_transient_correlation(da_act,this_ach,1,...
                'input_idx',eta_idx,'ref_idx_of_int',ref_idx_of_int,...
                'corr_idx_of_int',corr_idx,'local_transient_window',9);                                     
            null_r(i,:) = max(abs(this_corr.corr.corr_r));
        end
        for r = 1:size(da_act,2)
            this_null_r = transpose(null_r(:,r));
            this_corr = abs(mouse_corr.corr.corr_r(:,r));
            mouse_corr.corr.null_p(:,r) = sum(repmat(this_null_r,size(this_corr,1),1) > this_corr,2)/numel(this_null_r);        
        end
        rew_da.(mouse) = mouse_corr;
    end
end
        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_ACh_peak_rew_cc(
function rew_ach = get_ACh_peak_rew_cc(rew_data,null_it,varargin)

    %%%  parse optional inputs %%%
    ip = inputParser;
    ip.addParameter('eta_idx',-18:27);
    ip.addParameter('ref_idx_of_int',-9:9);
    ip.addParameter('corr_idx',0:9); 
    ip.parse(varargin{:});
    for j=fields(ip.Results)'
        eval([j{1} '=ip.Results.' j{1} ';']);
    end
    
    % loop
    mice = fieldnames(rew_data);    
    for m = 1:numel(mice)
        mouse = mice{m};   
        disp(mouse)
    %     mouse_idx = ismember(fib.mouse,mouse);
        da_act = rew_data.(mouse).DA.activity;
        ach_act = rew_data.(mouse).ACh.activity;
        % correlate ACh with DA peak magnitude
        mouse_corr = raster_transient_correlation(ach_act,da_act,1,...
            'input_idx',eta_idx,'ref_idx_of_int',ref_idx_of_int,...
            'corr_idx_of_int',corr_idx,'local_transient_window',9); 
        null_r = nan(null_it,size(da_act,2));
        % now null: shuffle the trials
        for i = 1:null_it      
            if rem(i,500) == 0
                disp(['     ' num2str(i)])
            end
            this_da = da_act(:,:,randperm(size(da_act,3)));
            this_corr = raster_transient_correlation(ach_act,this_da,1,...
                'input_idx',eta_idx,'ref_idx_of_int',ref_idx_of_int,...
            'corr_idx_of_int',corr_idx,'local_transient_window',9); 
            null_r(i,:) = max(abs(this_corr.corr.corr_r));
        end
        for r = 1:size(da_act,2)
            this_null_r = transpose(null_r(:,r));
            this_corr = abs(mouse_corr.corr.corr_r(:,r));
            mouse_corr.corr.null_p(:,r) = sum(repmat(this_null_r,size(this_corr,1),1) > this_corr,2)/numel(this_null_r);        
        end
        rew_ach.(mouse) = mouse_corr;
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
     
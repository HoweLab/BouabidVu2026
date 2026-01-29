%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organization
addpath(fullfile(pwd,'common_functions'))
data_dir = 'D:';
mice = {'UG27','UG28','UG29','UG30','UG31'};
fib = cohort_fib_table(data_dir,mice);
save_dir1 = fullfile(data_dir,'results','1_cross_corr');
corr_hotspot = load(fullfile(save_dir1,'cross_corr_dominant_results.mat'),'sig_moran','str');
save_dir3 = fullfile(data_dir,'results','3_unpred_rew');
if ~exist(save_dir3,'dir')
    mkdir(save_dir3)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. get concatenated reward data 
rew_data = get_all_rew_data(mice,fib,data_dir); % useful to save this


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. cross correlation of trial-by-trial ACh w trial-by-trial DA peak 
rew_cc = get_all_rew_cc(rew_data,fib);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. maps

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
results.vals.r = rew_cc.dominant.corr_r;
results.vals.p = rew_cc.dominant.p;
results.vals.lat = rew_cc.dominant.lat;
results.vals.dff_ach = rew_cc.dominant.dff_ach;
results.vals.cluster = vec(kmeans([results.vals.r,results.vals.lat],2));

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
save(fullfile(save_dir3,'unpred_rew_dominant_results.mat'),'-struct','results')
    
%%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. results figs
results = load(fullfile(save_dir3,'unpred_rew_dominant_results.mat'));

% scatter plot
plot_lag_corr_scatter(results.vals.r,results.vals.lat/18*1000,...
    'lat_bins',[-500:(1000/18):500],'clust_id',results.vals.cluster);

% smooth maps
ur_outlines = get_mask_projection_outlines(results.sig_moran.map,...
    results.str,'apply_str_mask',1,'proj_orientations',{'axial','sagittal'});
corr_outlines = get_mask_projection_outlines(corr_hotspot.sig_moran.map,...
    corr_hotspot.str,'apply_str_mask',1,'proj_orientations',{'axial','sagittal'});
outlines.axial = [ur_outlines.axial; corr_outlines.axial];
outlines.sagittal = [ur_outlines.sagittal; corr_outlines.sagittal];
plot_smooth_maps(results.smooth,results.str,'outlines',outlines);

% in-v-out violin
plot_violin_in_out(results.vals.r,fib,corr_hotspot.str,corr_hotspot.sig_moran.map)

% venn diagrams of hotspot comparisons (title has #voxels)
plot_hotspot_comparison_venn(results.sig_moran.corr_hotspot_overlap)
    

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_all_rew_data
function rew_data = get_all_rew_data(mice,fib,data_dir)
    neuromods = {'DA','ACh'};
    eta_idx = -9:27; % -.5s to 1.5s
    rew_data = struct;
    for m = 1:numel(mice)
        mouse = mice{m};           
        str_rois = fib.ROI_orig(strcmp(fib.mouse,mouse));

        % loop over data, keeping only unpred rew sessions, i.e.,
        % folders end in 'r'
        exp_dirs = dir(fullfile(data_dir,mouse));
        is_dirs = [exp_dirs.isdir];
        exp_dirs = {exp_dirs.name}';
        exp_dirs = exp_dirs(is_dirs);
        exp_dirs = exp_dirs(~startsWith(exp_dirs,'.') & endsWith(exp_dirs,'r'));
        
        % tmp cells for reward triggered averages
        tmp = struct;
        tmp.DA = cell(0,0);
        tmp.ACh = cell(0,0);   
        for d = 1:numel(exp_dirs)
            data = load(fullfile(data_dir,mouse,exp_dirs{d},[mouse '_' exp_dirs{d} '.mat']));
            if ~isempty(data.ACh.Fc_exp_hp_art) && ~isempty(data.DA.Fc_exp_hp_art)
                for n = 1:numel(neuromods)
                    rew_info = get_reward_info(data.(['behav_' neuromods{n}]));
                    these_rew = rew_info.rew_lick(~isnan(rew_info.rew_lick)); % trigger on consumption
                    eta = eventTriggeredAverage(data.(neuromods{n}).Fc_exp_hp_art(:,str_rois), ...
                        these_rew, eta_idx(1),eta_idx(end),'nullDistr',1,'bootstrapN',1);
                    tmp.(neuromods{n})(end+1,:) = [{eta} {data.(neuromods{n}).Fc_exp_hp_art(:,str_rois)}];
                end
            end   
        end
        for n = 1:numel(neuromods)
            eta = eta_combine(tmp.(neuromods{n}));
            eta = eta_sig(eta,'sig_idx',[find(eta_idx==0) find(eta_idx==18)]); % significance assessed 0-1s
            rew_data.(mouse).(neuromods{n}) = eta;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_all_cc
function rew_cc = get_all_rew_cc(rew_data,fib)
mice = fieldnames(rew_data);
eta_idx = -9:27; % -.5s to 1.5s   
which_r = {'min','max','dominant'};

% initialize
for i = 1:numel(which_r)
    rew_cc.(which_r{i}).corr_r = nan(size(fib,1),1);
    rew_cc.(which_r{i}).p = nan(size(fib,1),1);
    rew_cc.(which_r{i}).lat = nan(size(fib,1),1);
    rew_cc.(which_r{i}).dff_ach = nan(size(fib,1),1);
end

for m = 1:numel(mice)
    mouse = mice{m};   
    mouse_idx = ismember(fib.mouse,mouse);
    da_act = rew_data.(mouse).DA.activity;
    ach_act = rew_data.(mouse).ACh.activity;
    % correlate ACh with DA peak magnitude
    output = raster_transient_correlation(da_act,ach_act,1,...
        'input_idx',eta_idx,'ref_idx_of_int',0:18,'corr_idx_of_int',-9:9,...
        'local_transient_window',9,'center_trial_transient',1); 
    for i = 1:numel(which_r)
        rew_cc.(which_r{i}).corr_r(mouse_idx) = output.corr.([which_r{i} '_corr_r']);
        rew_cc.(which_r{i}).p(mouse_idx) = output.corr.([which_r{i} '_p']);
        rew_cc.(which_r{i}).lat(mouse_idx) = output.corr.([which_r{i} '_corr_r_input_idx']);       
        rew_cc.(which_r{i}).dff_ach(mouse_idx) = output.corr_raster.mu(...
            sub2ind(size(output.corr_raster.mu),...
            output.corr.([which_r{i} '_corr_r_idx']),...
            vec(1:size(output.corr_raster.mu,2))));
    end
end
end
        
     
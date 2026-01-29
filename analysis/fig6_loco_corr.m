%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organization
addpath(fullfile(pwd,'common_functions'))
data_dir = 'D:';
mice = {'UG27','UG28','UG29','UG30','UG31'};
fib = cohort_fib_table(data_dir,mice);
% load corr hotspot
save_dir1 = fullfile(data_dir,'results','1_cross_corr');
corr_hotspot = load(fullfile(save_dir1,'cross_corr_dominant_results.mat'),'sig_moran','str');
% load loco dynamics results
save_dir5 = fullfile(data_dir,'results','5_loco_dynamics');
trg_avgs = load(fullfile(save_dir5,'loco_trg_avg.mat'));
% directory for saving (interim) results
save_dir6 = fullfile(data_dir,'results','6_loco_corr');
if ~exist(save_dir6,'dir')
    mkdir(save_dir6)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. triggered average cross correlation
sr = 18; % sampling rate
neuromods = {'ACh','DA'};
loco_evs = {'init','invig'};
eta_idx = -18:36;
idx_of_int = -9:18;
loco_cc = get_all_loco_cc(trg_avgs,fib,'eta_idx',eta_idx,'idx_of_int',idx_of_int);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. maps
% size & striatum mask
str = trg_avgs.str;

for e = 1:numel(loco_evs)
    
    % sig_DA
    sig_da = abs(transpose(struct2array(structfun(@(x) x.DA_sig',...
        trg_avgs.(loco_evs{e}),'UniformOutput',false))));
    if strcmp(loco_evs{e},'init')
        sig_da = sig_da(:,2); % significant dip
    else
        sig_da = sig_da(:,1); % significant peak
    end
    
    % smooth maps and moran calculations
    results = struct;

    % vals
    results.vals.r = loco_cc.(loco_evs{e}).dominant.corr_r;
    results.vals.p = loco_cc.(loco_evs{e}).dominant.p;
    results.vals.lat = loco_cc.(loco_evs{e}).dominant.lat;
    results.vals.dff_ach = loco_cc.(loco_evs{e}).dominant.dff_ach;
    results.vals.sig_da = sig_da;
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
    save(fullfile(save_dir6,[loco_evs{e} '_dominant_results.mat']),'-struct','results')
end

%%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. results figs
for e = 1:numel(loco_evs)
    results = load(fullfile(save_dir6,[loco_evs{e} '_dominant_results.mat']));
    this_fib = fib(results.vals.sig_da==1,:);
    this_fib.ROI = vec(1:size(this_fib,1));
    
    % scatter plot
    plot_lag_corr_scatter(results.vals.r(results.vals.sig_da==1),...
        results.vals.lat(results.vals.sig_da==1)/18*1000,...
        'lat_bins',[-500:(1000/18):500],'clust_id',...
        results.vals.cluster(results.vals.sig_da==1));

    % smooth maps
    loco_outlines = get_mask_projection_outlines(results.sig_moran.map,...
        results.str,'apply_str_mask',1,'proj_orientations',{'axial','sagittal'});
    corr_outlines = get_mask_projection_outlines(corr_hotspot.sig_moran.map,...
        corr_hotspot.str,'apply_str_mask',1,'proj_orientations',{'axial','sagittal'});
    outlines.axial = [loco_outlines.axial; corr_outlines.axial];
    outlines.sagittal = [loco_outlines.sagittal; corr_outlines.sagittal];
    plot_smooth_maps(results.smooth,results.str,'outlines',outlines);

    % in-v-out violin
    plot_violin_in_out(results.vals.r(results.vals.sig_da==1),...
        this_fib,corr_hotspot.str,corr_hotspot.sig_moran.map)
    
    % venn diagrams of hotspot comparisons (title has #voxels)
    plot_hotspot_comparison_venn(results.sig_moran.corr_hotspot_overlap)

end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_all_loco_cc
function loco_cc = get_all_loco_cc(trg_avgs,fib,varargin)

%%%  parse optional inputs %%%
ip = inputParser;
ip.addParameter('eta_idx',-18:36)
ip.addParameter('idx_of_int',-9:18)
ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end

% initialize
mice = unique(fib.mouse);
which_r = {'min','max','dominant'};
loco_evs = {'init','invig'};
loco_cc = struct;
for e = 1:numel(loco_evs)
    for i = 1:numel(which_r)
        loco_cc.(loco_evs{e}).(which_r{i}).corr_r = nan(size(fib,1),1);
        loco_cc.(loco_evs{e}).(which_r{i}).p = nan(size(fib,1),1);
        loco_cc.(loco_evs{e}).(which_r{i}).lat = nan(size(fib,1),1);
        loco_cc.(loco_evs{e}).(which_r{i}).dff_ach = nan(size(fib,1),1);
    end
end

for e = 1:numel(loco_evs)
    if strcmp(loco_evs{e},'init')
        DA_tr_sign = -1; % dip
    else
        DA_tr_sign = 1; % peak
    end
    for m = 1:numel(mice)
        mouse = mice{m};
        mouse_idx = ismember(fib.mouse,mouse);
        da_act = trg_avgs.(loco_evs{e}).(mouse).DA;
        ach_act = trg_avgs.(loco_evs{e}).(mouse).ACh;
        output = raster_transient_correlation(da_act,ach_act,DA_tr_sign,...
            'input_idx',eta_idx,'ref_idx_of_int',idx_of_int,'corr_idx_of_int',-9:9,...
            'local_transient_window',9,'center_trial_transient',1); 
 
        for i = 1:numel(which_r)
            loco_cc.(loco_evs{e}).(which_r{i}).corr_r(mouse_idx) = output.corr.([which_r{i} '_corr_r']);
            loco_cc.(loco_evs{e}).(which_r{i}).p(mouse_idx) = output.corr.([which_r{i} '_p']);
            loco_cc.(loco_evs{e}).(which_r{i}).lat(mouse_idx) = output.corr.([which_r{i} '_corr_r_input_idx']);       
            loco_cc.(loco_evs{e}).(which_r{i}).dff_ach(mouse_idx) = ...
                output.corr_raster.mu(...
                sub2ind(size(output.corr_raster.mu),...
                output.corr.([which_r{i} '_corr_r_idx']),...
                vec(1:size(output.corr_raster.mu,2))));
        end
    end
end
end
    

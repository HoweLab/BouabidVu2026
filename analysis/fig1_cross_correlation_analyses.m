%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organization
addpath(fullfile(pwd,'common_functions'))
data_dir = 'D:';
mice = {'UG27','UG28','UG29','UG30','UG31'};
fib = cohort_fib_table(data_dir,mice);

% directory for saving (interim) results
save_dir1 = fullfile(data_dir,'results','1_cross_corr');
if ~exist(save_dir1,'dir')
    mkdir(save_dir1)
end
sr = 18; % sampling rate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. concatenate data 
nanpad = sr*5; % put a bunch of nans between days so the cross-correlation sliding doesn't combine data across days
concat_data = get_concat_data(mice,fib,data_dir,nanpad);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. cross-correlation
% note: lags = ACh lag w/respect to DA, i.e. (+)lag means DA leads
lag = 18;
cc_results = get_cross_corr_results(concat_data,mice,fib,lag);
clear concat_data % clear for memory; no longer need

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

% save for convenience
results.str = str;
save(fullfile(save_dir1,'cross_corr_dominant_results.mat'),'-struct','results')
    
%%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. results figs
results = load(fullfile(save_dir1,'cross_corr_dominant_results.mat'));

% scatter plot
plot_lag_corr_scatter(results.vals.r,results.vals.lag/18*1000,'lat_bins',[-1000:(1000/9):1000]);

% smooth maps
outlines = get_mask_projection_outlines(results.sig_moran.map,...
    results.str,'apply_str_mask',1,'proj_orientations',{'axial','sagittal'});    
plot_smooth_maps(results.smooth,results.str,'outlines',outlines);

% in-v-out violin
plot_violin_in_out(results.vals.r,fib,results.str,results.sig_moran.map)

% pie chart: mouse composition of hotspot
map_ind = get_map_ind(fib,results.str);
in_map = results.sig_moran.map(map_ind) == 1;
mouse_n = zeros(numel(mice),1);
hotspot_fib = fib(in_map,:);
for m = 1:numel(mice)
    mouse_n(m) =sum(ismember(cellstr(hotspot_fib.mouse),mice{m}));
end
figure
p = pie(mouse_n);
p_colors = lines(7);
p_colors(4,:) = p_colors(6,:);
for i = 1:numel(mice)
    p((2*i)-1).FaceColor = p_colors(i,:);
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_concat_data
function concat_data = get_concat_data(mice,fib,data_dir,nanpad)
    concat_data = struct;
    for m = 1:numel(mice)
        mouse = mice{m};           
        str_rois = fib.ROI_orig(strcmp(fib.mouse,mouse));

        % preallocate
        concat_data.(mouse).date = [];
        concat_data.(mouse).DA = [];
        concat_data.(mouse).ACh = [];

        % loop over data, keeping only spontaneous movement sessions, i.e.,
        % non-reward sessions (non-reward session folders end in 'r')
        exp_dirs = dir(fullfile(data_dir,mouse));
        is_dirs = [exp_dirs.isdir];
        exp_dirs = {exp_dirs.name}';
        exp_dirs = exp_dirs(is_dirs);
        exp_dirs = exp_dirs(~startsWith(exp_dirs,'.') & ~endsWith(exp_dirs,'r'));
        for d = 1:numel(exp_dirs)
            data = load(fullfile(data_dir,mouse,exp_dirs{d},[mouse '_' exp_dirs{d} '.mat']));
            if ~isempty(data.ACh.Fc_exp_hp_art) && ~isempty(data.DA.Fc_exp_hp_art)
                % append
                concat_data.(mouse).DA = [concat_data.(mouse).DA; ...
                    nan(nanpad,numel(str_rois)); ...
                    data.DA.Fc_exp_hp_art(:,str_rois)];
                concat_data.(mouse).ACh = [concat_data.(mouse).ACh;...
                    nan(nanpad,numel(str_rois)); ...
                    data.ACh.Fc_exp_hp_art(:,str_rois)];
                concat_data.(mouse).date = [concat_data.(mouse).date;...
                    repmat({'n/a'},nanpad,1); ...
                    repmat(exp_dirs(d),size(data.DA.Fc_exp_hp_art,1),1)];
            end   
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_cross_corr_results
function cc_results = get_cross_corr_results(concat_data,mice,fib,lag)
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

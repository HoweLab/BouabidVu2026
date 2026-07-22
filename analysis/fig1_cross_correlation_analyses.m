%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organization
addpath(fullfile(pwd,'common_functions'))
data_dir = 'G:';
mice = {'UG27','UG28','UG29','UG30','UG31'};
fib = cohort_fib_table(data_dir,mice);

% directory for saving (interim) results
save_dir1 = fullfile(data_dir,'results','1_cross_corr');
if ~exist(save_dir1,'dir')
    mkdir(save_dir1)
end
sr = 18; % sampling rate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. cross-correlation across sessions & relevant lags
lag = sr; % lag = +/- one second (sampling rate)

% 1a. cross-correlation across sessions
[cc_results,mouse_results] = get_cross_corr_results(mice,fib,data_dir,lag,...
   'save_by_mouse',1,'save_dir',save_dir1,'null_n',10000);

% 1b. distribution of relevant lags and cross correlation at those lags
cc_lags_distr = get_cross_corr_lag_distr(mice,fib,data_dir,lag,mouse_results,...
    'save_dir',save_dir1,'sample_n',10000);

% 1c. histogram of lags and identification of components
f1 = plot_lag_hist_and_components(cc_lags_distr,lag);

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. maps: interpolated

% load if available
map_results = load_if_exist(fullfile(save_dir1,'map_results.mat'));

% size & striatum mask
if ~isfield(map_results,'str')
    voxel_size = 0.05;
    str = get_striatum_vol_mask(...
        [min(fib.fiber_bottom_AP) max(fib.fiber_bottom_AP)],... % AP_range
        [min(fib.fiber_bottom_ML) max(fib.fiber_bottom_ML)],... % ML_range
        [min(fib.fiber_bottom_DV) max(fib.fiber_bottom_DV)],... % DV_range
        voxel_size);
    map_results.str = str;
end

% vals
map_results.vals.r.neg_lag = tanh(cc_lags_distr.lag_gm.all.weighted_r_z(:,cc_lags_distr.lag_gm.all.main_neg_idx));
map_results.vals.r.pos_lag = tanh(cc_lags_distr.lag_gm.all.weighted_r_z(:,cc_lags_distr.lag_gm.all.main_pos_idx));

% smooth maps
if ~isfield(map_results,'info')
    lag_signs = {'neg','pos'};
    tmp = get_activity_map_interp([map_results.vals.r.neg_lag map_results.vals.r.pos_lag],...
        fib, str.info.voxel_size,'AP_range',[min(str.info.AP) max(str.info.AP)],...
        'ML_range',[min(str.info.ML) max(str.info.ML)],...
        'DV_range',[min(str.info.DV) max(str.info.DV)],'incl_plot_info',1);  
    for i = 1:numel(lag_signs)
        map_results.interp.([lag_signs{i} '_lag']).vol = tmp.(['vol_' sprintf('%02d',i)]).interp;   % interpolated volume
        map_results.interp.([lag_signs{i} '_lag']).n = tmp.(['vol_' sprintf('%02d',i)]).n_mice;     % #mice contrib to each voxel
        map_results.interp.([lag_signs{i} '_lag']).F = tmp.(['vol_' sprintf('%02d',i)]).interp_F;   % interpolant function
    end
    map_results.info = tmp.info;
    save(fullfile(save_dir1,'map_results.mat'),'-struct','map_results','-v7.3')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. maps: moran
lag_signs = {'neg','pos'};

% 3a. first determine proper neighborhood cube width
% 3ai. first figure out a distance based on each mouse's coverage
if ~isfield(map_results,'moran') || ~isfield(map_results.moran,'neighborhood') || ...
        ~isfield(map_results.moran.neighborhood,'mouse_neighborhood')
    cohort_min_fib = round(mean(cellfun(@(x) sum(ismember(fib.mouse,x))^(1/3),mice)));
    tmp = get_mouse_neighborhood_r(fib,'min_n_fibs',cohort_min_fib);
    map_results.moran.neighborhood.mouse_neighborhood.n_fibs = cohort_min_fib;
    map_results.moran.neighborhood.mouse_neighborhood.mm = max(structfun(@(x) x.r,tmp));
    map_results.moran.neighborhood.mouse_neighborhood.vox = round(max(structfun(@(x) x.r,tmp))/map_results.str.info.voxel_size);
    save(fullfile(save_dir1,'map_results.mat'),'-struct','map_results','-v7.3')
end

% 3aii. now based on a quick scan of local moran's
% let's test the values around the mouse neighborhood and go up to 31
widths_to_test = 3:2:31;
if ~isfield(map_results,'moran') || ~isfield(map_results.moran,'neighborhood') || ...
    ~isfield(map_results.moran.neighborhood,'moran_neighborhood')
    tmp = struct;
    for i = 1:numel(lag_signs)
        value_array = map_results.vals.r.([lag_signs{i} '_lag']);    
        tmp.([lag_signs{i} '_lag']) = get_moran_neighborhood_width(...
            fib,value_array,map_results.str.info.voxel_size,widths_to_test,'n_it',500,...
            'AP_range',[min(map_results.str.info.AP) max(map_results.str.info.AP)],...
            'ML_range',[min(map_results.str.info.ML) max(map_results.str.info.ML)],...
            'DV_range',[min(map_results.str.info.DV) max(map_results.str.info.DV)]);
        tmp.([lag_signs{i} '_lag']).loc_max = tmp.([lag_signs{i} '_lag']).width(...
            find(islocalmax(tmp.([lag_signs{i} '_lag']).width_z),1,'first'));
    end
    map_results.moran.neighborhood.moran_neighborhood = tmp;
    save(fullfile(save_dir1,'map_results.mat'),'-struct','map_results','-v7.3')
end

% 3b. actual moran, based on that neighborhood
for i = 1:numel(lag_signs)
    if ~isfield(map_results.moran,[lag_signs{i} '_lag'])
        disp([lag_signs{i} '_lag: calculating local Moran''s I and null distr']);
        neighborhood_width = max([...
            map_results.moran.neighborhood.moran_neighborhood.([lag_signs{i} '_lag']).loc_max,...
            map_results.moran.neighborhood.mouse_neighborhood.vox,...
            ]);
        map_results.moran.neighborhood.width.([lag_signs{i} '_lag']) = neighborhood_width; % store this
        data = map_results.interp.([lag_signs{i} '_lag']).vol;
        map_results.moran.([lag_signs{i} '_lag']) = ...
            local_morans_I_bootstrap_null(...
            fib,map_results.vals.r.([lag_signs{i} '_lag']),...
            map_results.str.info.voxel_size,'n_it',10000,'batch_size',1000,...
            'AP_range',[min(map_results.str.info.AP) max(map_results.str.info.AP)],...
            'ML_range',[min(map_results.str.info.ML) max(map_results.str.info.ML)],...
            'DV_range',[min(map_results.str.info.DV) max(map_results.str.info.DV)],...
            'prctiles',[95 97.5 99 99.5 99.9],'save_null',fullfile(save_dir1,['null_moran_' lag_signs{i} '_lag.mat']),...
            'weight_matrix',ones(neighborhood_width,neighborhood_width,neighborhood_width));
        save(fullfile(save_dir1,'map_results.mat'),'-struct','map_results','-v7.3')
    end
end
% 3c. significant hotspot
for i = 1:numel(lag_signs)
    if ~isfield(map_results.moran.([lag_signs{i} '_lag']),'sig')    
        map_results.moran.([lag_signs{i} '_lag']) = local_morans_I_sig_moran(...
            map_results.moran.([lag_signs{i} '_lag']),...
            fullfile(save_dir1,['null_moran_' lag_signs{i} '_lag.mat']),...
            'mask',str.striatum_mask);  
        save(fullfile(save_dir1,'map_results.mat'),'-struct','map_results','-v7.3')
    end
end

    
% %%   
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % 4. results figs
% % results = load(fullfile(save_dir1,'cross_corr_dominant_results.mat'));
% % 
% % % scatter plot
% % plot_lag_corr_scatter(results.vals.r,results.vals.lag/18*1000,'lat_bins',[-1000:(1000/9):1000]);
% % 
% % % smooth maps
% % outlines = get_mask_projection_outlines(results.sig_moran.map,...
% %     results.str,'apply_str_mask',1,'proj_orientations',{'axial','sagittal'});    
% % plot_smooth_maps(results.smooth,results.str,'outlines',outlines);
% % 
% % % in-v-out violin
% % plot_violin_in_out(results.vals.r,fib,results.str,results.sig_moran.map)
% % 
% % % pie chart: mouse composition of hotspot
% % map_ind = get_map_ind(fib,results.str);
% % in_map = results.sig_moran.map(map_ind) == 1;
% % mouse_n = zeros(numel(mice),1);
% % hotspot_fib = fib(in_map,:);
% % for m = 1:numel(mice)
% %     mouse_n(m) =sum(ismember(cellstr(hotspot_fib.mouse),mice{m}));
% % end
% % figure
% % p = pie(mouse_n);
% % p_colors = lines(7);
% % p_colors(4,:) = p_colors(6,:);
% % for i = 1:numel(mice)
% %     p((2*i)-1).FaceColor = p_colors(i,:);
% % end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_cross_corr_results
function [cc_results,mouse_results] = get_cross_corr_results(mice,fib,data_dir,lag,varargin)
    
    %%%  parse optional inputs %%%
    ip = inputParser;
    ip.addParameter('null_n',10000);        % # iterations to run for each session's null distribution
    ip.addParameter('alpha_val',5);         % alpha value in %
    ip.addParameter('save_by_mouse',1);     % save out each mouse's data
    ip.addParameter('save_dir',[]);         % save directory
    ip.parse(varargin{:});
    for j=fields(ip.Results)'
        eval([j{1} '=ip.Results.' j{1} ';']);
    end
    
    % if there's no save dir we don't save anything
    if isempty(save_dir)
        save_by_mouse = 0;
    end

    if ~isempty(save_dir) && exist(fullfile(save_dir,'cross_corr_results.mat'),'file')
        cc_results = load(fullfile(save_dir,'cross_corr_results.mat'));
        mouse_results = load(fullfile(save_dir,'cross_corr_mouse_results.mat'));        
    else
        % initialize output
        cc_results = struct;                
        cc_results.lag = -lag:lag;
        cc_results.mouse = [];
        cc_results.str_rois = [];
        cc_results.r.mean = [];
        cc_results.r.std = [];
        cc_results.r.n = [];
        cc_results.r.sig = [];
        cc_results.r.max = [];
        cc_results.r.min = [];
        mouse_results = struct;

        % get per session cross-correlation and null cross-correlation 
        for m = 1:numel(mice)        
            mouse = mice{m};
            disp(mouse)        
            if save_by_mouse == 1 && exist(fullfile(save_dir,[mouse '.mat']),'file')
                session_cross_corrs = load(fullfile(save_dir,[mouse '.mat']));
            else
                session_cross_corrs = get_session_cross_corrs(mouse,fib,data_dir,lag);
                session_cross_corrs.null = get_null_cross_corrs(mouse,fib,data_dir,lag,'n_it',null_n);
                session_cross_corrs.null = session_cross_corrs.null.max_abs_r; 
                if save_by_mouse == 1
                    save(fullfile(save_dir,[mouse '.mat']),'-struct','session_cross_corrs','-v7.3')
                end
            end
            mouse_results.(mouse) = session_cross_corrs;
            
            % extract useful statistics and add to results
            cc_stats = get_session_cross_corr_stats(session_cross_corrs,alpha_val);
            cc_results.mouse = [cc_results.mouse; repmat({mouse},numel(session_cross_corrs.str_rois),1)];
            cc_results.str_rois = [cc_results.str_rois; vec(session_cross_corrs.str_rois)];
            cc_stats_fields = fieldnames(cc_results.r);
            for f = 1:numel(cc_stats_fields)
                cc_results.r.(cc_stats_fields{f}) = ...
                    [cc_results.r.(cc_stats_fields{f});cc_stats.(cc_stats_fields{f})];
            end
        end

        % save
        if ~isempty(save_dir)
            save(fullfile(save_dir,'cross_corr_results.mat'),'-struct','cc_results')
            save(fullfile(save_dir,'cross_corr_mouse_results.mat'),'-struct','mouse_results','-v7.3')
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_cross_corr_lag_distr:
function cc_lags_distr = get_cross_corr_lag_distr(mice,fib,data_dir,lag,...
    mouse_results,varargin)
    
    %%%  parse optional inputs %%%
    ip = inputParser;
    ip.addParameter('sample_n',10000);      % #iterations for bootstrapped distribution
    ip.addParameter('sample_n_chunk',500)   % #iterations to run at a time, for iterative saving
    ip.addParameter('n_calculation',[]);    % parameters for calculating min n needed for corr estimation    
    ip.addParameter('save_dir',[]);         % save directory
    ip.parse(varargin{:});
    for j=fields(ip.Results)'
        eval([j{1} '=ip.Results.' j{1} ';']);
    end
    
    % for the n_calculation
    if isempty(n_calculation)
        n_calculation = struct;
        n_calculation.r_target = 0.3;
        n_calculation.alpha = 0.05;
        n_calculation.power = 0.8;
    end
    
    % load if we have it
    cc_lags_distr = load_if_exist(fullfile(save_dir,'cross_corr_r_lag_results.mat'));      
    
    % each mouse's min #timepoints needed to estimate corr
    if ~isfield(cc_lags_distr,'min_n')
        % get per-mouse min #samples needed to detect weakish (0.3) corr
        for m = 1:numel(mice)        
            mouse = mice{m};
            disp(mouse)        
            % info
            exp_dirs = get_exp_dirs(mouse,data_dir);
            str_rois = fib.ROI_orig(strcmp(fib.mouse,mouse));
            
            % get min number of timepoints needed
            min_n = nan(numel(exp_dirs),numel(str_rois));
            for d = 1:numel(exp_dirs)
                exp_dir = exp_dirs{d};
                data = load(fullfile(data_dir,mouse,exp_dir,[mouse '_' exp_dir '.mat'])); % load data
                for roi = 1:numel(str_rois)
                    ACh = data.ACh.Fc_exp_hp_art(:,str_rois(roi));
                    DA = data.DA.Fc_exp_hp_art(:,str_rois(roi));
                    min_n(d,roi) = get_required_corr_n(DA,ACh,n_calculation.r_target,n_calculation.alpha,n_calculation.power,2*lag+1);
                end
            end
            min_n = ceil(max(min_n(:))); % most conservative estimate of min n needed;
            disp(['     min n: ' num2str(min_n)])
            cc_lags_distr.min_n.(mouse) = min_n;
        end        
    end
    
    % initialize output
    output_fields = {'gm_mu','gm_std','gm_perc','gm_r','gm_null_r','gm_r_sig',...
        'mouse','roi_num','loc_extr_sign','main_comp'};    
    if ~isfield(cc_lags_distr,'lag_hist')
        cc_lags_distr.lag_hist = struct;
    end
    
    % if we need to compile the data
    if ~isfield(cc_lags_distr.lag_hist,output_fields{end})
        % initialize
        for f = 1:numel(output_fields)                
            cc_lags_distr.lag_hist.(output_fields{f}) = [];
        end            
        % get per session cross-correlation and null cross-correlation 
        loc_signs = {'min','max'};
        for m = 1:numel(mice)        
            mouse = mice{m};
            disp(mouse)        
            % mouse lag data
            mouse_data = get_mouse_cross_corr_lag_distr(mouse,fib,data_dir,lag,cc_lags_distr.min_n.(mouse),'save_dir',save_dir);        
            % mouse cross corr data
            cc_mouse = mouse_results.(mouse);
            % add to results
            for s = 1:numel(loc_signs)
                for r = 1:numel(mouse_data.lags.(['loc_' loc_signs{s}]).gm_sig)
                    these_sig = find(mouse_data.lags.(['loc_' loc_signs{s}]).gm_sig{r});
                    disp(['     ' loc_signs{s} ' ' num2str(r)])
                    for i = 1:numel(these_sig)
                        gm_idx = these_sig(i);
                        % mean of this component
                        cc_lags_distr.lag_hist.gm_mu = [cc_lags_distr.lag_hist.gm_mu; ...
                            mouse_data.lags.(['loc_' loc_signs{s}]).gm{r}.mu(gm_idx)];
                        % std of this component
                        cc_lags_distr.lag_hist.gm_std = [cc_lags_distr.lag_hist.gm_std; ...
                            mouse_data.lags.(['loc_' loc_signs{s}]).gm{r}.Sigma(gm_idx)];
                        % percent of this component
                        cc_lags_distr.lag_hist.gm_perc = [cc_lags_distr.lag_hist.gm_perc; ...
                            mouse_data.lags.(['loc_' loc_signs{s}]).gm{r}.ComponentProportion(gm_idx)];                
                        % mouse
                        cc_lags_distr.lag_hist.mouse = [cc_lags_distr.lag_hist.mouse;...
                            {mouse}];
                        % str_roi num
                        cc_lags_distr.lag_hist.roi_num = [cc_lags_distr.lag_hist.roi_num;...
                            cc_mouse.str_rois(r)];
                        % local extrema sign
                        cc_lags_distr.lag_hist.loc_extr_sign = ...
                            [cc_lags_distr.lag_hist.loc_extr_sign; loc_signs{s}];
                        % whether it's the main comp (based on percent)
                        [~,is_main] = max(mouse_data.lags.(['loc_' loc_signs{s}]).gm{r}.ComponentProportion);
                        cc_lags_distr.lag_hist.main_comp = ...
                            [cc_lags_distr.lag_hist.main_comp; gm_idx==is_main];
                        % weighted r    
                        gm_r = get_gmm_weighted_mean(permute(cc_mouse.r(r,:,:),[3 2 1]),...
                            mouse_data.lags.(['loc_' loc_signs{s}]).gm{r},-lag:lag);                
                        cc_lags_distr.lag_hist.gm_r = [cc_lags_distr.lag_hist.gm_r; mean(gm_r(:,gm_idx))];                    
                        % weighted null r 95%
                        null_r_95 = repmat(prctile(cc_mouse.null(r,:,:),95,2),1,numel(-lag:lag),1);
                        gm_null_r = get_gmm_weighted_mean(permute(null_r_95,[3 2 1]),...
                            mouse_data.lags.(['loc_' loc_signs{s}]).gm{r},-lag:lag);                                  
                        cc_lags_distr.lag_hist.gm_null_r = [cc_lags_distr.lag_hist.gm_null_r;mean(gm_null_r(:,gm_idx))];
                        % significance
                        cc_lags_distr.lag_hist.gm_r_sig = [cc_lags_distr.lag_hist.gm_r_sig;...
                        abs(mean(gm_r(:,gm_idx))) > mean(gm_null_r(:,gm_idx))];
                    end
                end
            end                
        end
        if ~isempty(save_dir)
            save(fullfile(save_dir,'cross_corr_r_lag_results.mat'),'-struct','cc_lags_distr')
        end
    end
    
    % fit overall GM
    if ~isfield(cc_lags_distr,'lag_gm') || ~isfield(cc_lags_distr.lag_gm,'all')
        [this_gm,this_gof] = fit_gmm_to_hist(cc_lags_distr.lag_hist.gm_mu,...
            'gof','BIC','choose','elbow','max_n',18,'n_tries',5);
        cc_lags_distr.lag_gm.all.gm = this_gm;
        cc_lags_distr.lag_gm.all.gof = this_gof;
        cc_lags_distr.lag_gm.all.gm_r = get_gmm_weighted_mean(permute(null_r_95,[3 2 1]),...
            mouse_data.lags.(['loc_' loc_signs{s}]).gm{r},-lag:lag); 
        [cc_lags_distr.lag_gm.all.main_neg_idx, cc_lags_distr.lag_gm.all.main_pos_idx] = ...
            get_main_neg_pos_gmm_idx(cc_lags_distr.lag_gm.all.gm);
    end
    
    % get weighted means for each mouse
    cc_lags_distr.lag_gm.all.weighted_r = [];
    cc_lags_distr.lag_gm.all.weighted_null = [];
    cc_lags_distr.lag_gm.all.weighted_r_z = [];
    cc_lags_distr.lag_gm.all.weighted_null_z = [];
    for m = 1:numel(mice)
        mouse = mice{m};
        mouse_cc = mouse_results.(mouse);
        for r = 1:numel(mouse_cc.str_rois)
            this_weighted_r = get_gmm_weighted_mean(...
                permute(mouse_cc.r(r,:,:),[3 2 1]),cc_lags_distr.lag_gm.all.gm,-lag:lag);
            this_weighted_z = get_gmm_weighted_mean(...
                permute(atanh(mouse_cc.r(r,:,:)),[3 2 1]),cc_lags_distr.lag_gm.all.gm,-lag:lag);
            this_weighted_null = get_gmm_weighted_mean(...
                permute(repmat(prctile(mouse_cc.null(r,:,:),95,2),1,numel(-lag:lag),1),[3 2 1]),...
                cc_lags_distr.lag_gm.all.gm,-lag:lag);
            this_weighted_null_z = get_gmm_weighted_mean(...
                permute(repmat(prctile(atanh(mouse_cc.null(r,:,:)),95,2),1,numel(-lag:lag),1),[3 2 1]),...
                cc_lags_distr.lag_gm.all.gm,-lag:lag);
            
            cc_lags_distr.lag_gm.all.weighted_r = [cc_lags_distr.lag_gm.all.weighted_r;...
                mean(this_weighted_r)];
            cc_lags_distr.lag_gm.all.weighted_r_z = [cc_lags_distr.lag_gm.all.weighted_r_z;...
                mean(this_weighted_z)];
            cc_lags_distr.lag_gm.all.weighted_null = [cc_lags_distr.lag_gm.all.weighted_null;...
                mean(this_weighted_null)];
            cc_lags_distr.lag_gm.all.weighted_null_z = [cc_lags_distr.lag_gm.all.weighted_null_z;...
                mean(this_weighted_null_z)];
        end
    end
    cc_lags_distr.lag_gm.all.weighted_r_sig = ...
        abs(cc_lags_distr.lag_gm.all.weighted_r) > cc_lags_distr.lag_gm.all.weighted_null;
    cc_lags_distr.lag_gm.all.weighted_r_z_sig = ...
        abs(cc_lags_distr.lag_gm.all.weighted_r_z) > cc_lags_distr.lag_gm.all.weighted_null_z;
    
    % save
    if ~isempty(save_dir)
        save(fullfile(save_dir,'cross_corr_r_lag_results.mat'),'-struct','cc_lags_distr')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_mouse_cross_corr_lag_distr:
function mouse_data = get_mouse_cross_corr_lag_distr(mouse,fib,data_dir,lag,min_n,varargin)
    
    %%%  parse optional inputs %%%
    ip = inputParser;
    ip.addParameter('sample_n',10000);                  % #iterations for bootstrapped distribution
    ip.addParameter('sample_n_chunk',500)               % #iterations to run at a time, for iterative saving   
    ip.addParameter('n_calculation',[]);                % parameters for calculating min n needed for corr estimation
    ip.addParameter('plot_hists',0);                    % whether or not to plot the lag histograms (1 per site)
    ip.addParameter('save_dir',[]);                     % save directory
    ip.addParameter('save_suffix','_r_lag_distr.mat');  % save suffix
    ip.parse(varargin{:});
    for j=fields(ip.Results)'
        eval([j{1} '=ip.Results.' j{1} ';']);
    end
       
    % info
    exp_dirs = get_exp_dirs(mouse,data_dir);
    str_rois = fib.ROI_orig(strcmp(fib.mouse,mouse));

    % preallocate
    n = nan(sample_n,1);
    r = nan(numel(str_rois),2*lag+1,sample_n);
    max_null_r = nan(numel(str_rois),sample_n);
    
    % load if available; overwrite with existing data
    if ~isempty(save_dir) && exist(fullfile(save_dir,[mouse save_suffix]),'file')
        mouse_data = load(fullfile(save_dir,[mouse '_r_lag_distr.mat']));
        i_leftoff = find(~isnan(mouse_data.n),1,'last');
        n(1:i_leftoff) = mouse_data.n(1:i_leftoff);
        r(:,:,1:i_leftoff) = mouse_data.r(:,:,1:i_leftoff);
        max_null_r(:,1:i_leftoff) = mouse_data.max_null_r(:,1:i_leftoff);            
    end

    % for each iteration, randomly select a day and a chunk of time
    i_pickup = find(isnan(n),1,'first');
    if ~isempty(i_pickup) % if we don't have any to do
        while isnan(n(sample_n))
            disp(['     picking up from ' num2str(i_pickup)])
            % do this in chunks for easy iterative saving
            parfor i = i_pickup:min([i_pickup+sample_n_chunk sample_n])
                % randomly choose a date, time window, and time window start
                % taking into account the minumum n needed to estimate a pearson r of 0.3
                exp_dir = exp_dirs{randi([1 numel(exp_dirs)])}; % randomly choose a date                
                data = load(fullfile(data_dir,mouse,exp_dir,[mouse '_' exp_dir '.mat'])); % load data
                DA = data.DA.Fc_exp_hp_art(:,str_rois);
                ACh = data.ACh.Fc_exp_hp_art(:,str_rois);
                n_timepoints = randi([min_n size(DA,1)-2*lag]);
                i_start = randi([lag+1 size(DA,1)-n_timepoints-lag]);
                %disp(['          ' num2str([i i_start n_timepoints])])
                % get the cross correlation and null cross correlation 
                this_DA = DA(i_start:(i_start+n_timepoints-1),:);
                this_ACh = ACh(i_start:(i_start+n_timepoints-1),:);
                session_cross_corr = get_session_cross_corr(this_DA,this_ACh,lag);        
                null_r = get_session_cross_corr_null(DA,ACh,lag,'n_timepoints',n_timepoints,'n_it',1);                
                n(i) = n_timepoints;
                r(:,:,i) = session_cross_corr.r;
                max_null_r(:,i) = max(abs(null_r),[],2);
            end
            if ~isempty(save_dir)
                save(fullfile(save_dir,[mouse save_suffix]),'n','r','max_null_r','-v7.3')
            end
            i_pickup = i_pickup+sample_n_chunk;
        end                         
        mouse_data.n = n;
        mouse_data.r = r;
        mouse_data.max_null_r = max_null_r;                    
    end
    clear n r max_null_r

    % calculate p-value (if haven't)
    if ~isfield(mouse_data,'p')
        disp('     calculating p')
        mouse_data.p = nan(size(mouse_data.r));
        for roi_num = 1:size(mouse_data.r, 1)
            r_site   = abs(mouse_data.r(roi_num,:,:));                          % (1 x n_lags x n_sessions)
            null_site = permute(mouse_data.max_null_r(roi_num,:), [1 3 4 2]);   % (1 x 1 x 1 x n_perms)
            mouse_data.p(roi_num,:,:) = sum(null_site >= r_site, 4) / size(mouse_data.max_null_r,2);
        end            
    end  

    % get max and min (if haven't)
    if ~isfield(mouse_data,'lags') || ~isfield(mouse_data,'minmax')
        disp('     getting min & max r and lags')
        output =  get_cc_lag_sample_stats(mouse_data,0.05,-lag:lag);
        mouse_data.minmax.max = output.max;
        mouse_data.minmax.min = output.min;            
        mouse_data.lags.loc_max = output.loc_max;
        mouse_data.lags.loc_min = output.loc_min;                                                
    end        

    % now fit GMM to lags for pos corr and neg corr (if haven't)
    if ~isfield(mouse_data.lags,'loc_max') || ~isfield(mouse_data.lags,'loc_min')                           
        disp('     fitting gmm to lags')            
        loc_signs = {'min','max'};
        for s = 1:numel(loc_signs)                               
            mouse_data.lags.(['loc_' loc_signs{s}]).gm = cell(numel(str_rois),1);
            mouse_data.lags.(['loc_' loc_signs{s}]).gm_gof = cell(numel(str_rois),1);                
            mouse_data.lags.(['loc_' loc_signs{s}]).gm_weighted_r = cell(numel(str_rois),1);
            mouse_data.lags.(['loc_' loc_signs{s}]).gm_weighted_null_r = cell(numel(str_rois),1);
            mouse_data.lags.(['loc_' loc_signs{s}]).gm_sig = cell(numel(str_rois),1);
            for r = 1:numel(str_rois)    
                disp(['          ' loc_signs{s} ': ' num2str(r) ' of ' num2str(numel(str_rois))])
                X = mouse_data.lags.(['loc_' loc_signs{s}]).lags{r};                    
                [this_gm,this_gof] = fit_gmm_to_hist(X,'gof','BIC','choose','elbow','max_n',18,'n_tries',5);
                mouse_data.lags.(['loc_' loc_signs{s}]).gm{r} = this_gm;
                mouse_data.lags.(['loc_' loc_signs{s}]).gm_gof{r} = this_gof;
                % GMM-weighted r and null
                mouse_data.lags.(['loc_' loc_signs{s}]).gm_weighted_r{r} = ...
                    get_gmm_weighted_mean(permute(mouse_data.r(r,:,:),[3 2 1]),this_gm,-lag:lag);
                mouse_data.lags.(['loc_' loc_signs{s}]).gm_weighted_null_r{r} = ...
                    get_gmm_weighted_mean(permute(repmat(mouse_data.max_null_r(r,:,:),1,1,numel(-lag:lag)),[2 3 1]),this_gm,-lag:lag);
                if strcmp(loc_signs{s},'min')
                    mouse_data.lags.(['loc_' loc_signs{s}]).gm_sig{r} = ...
                        mean(mouse_data.lags.(['loc_' loc_signs{s}]).gm_weighted_r{r}) <...
                        -prctile(mouse_data.lags.(['loc_' loc_signs{s}]).gm_weighted_null_r{r},95);
                elseif strcmp(loc_signs{s},'max')
                    mouse_data.lags.(['loc_' loc_signs{s}]).gm_sig{r} = ...
                        mean(mouse_data.lags.(['loc_' loc_signs{s}]).gm_weighted_r{r}) >...
                        prctile(mouse_data.lags.(['loc_' loc_signs{s}]).gm_weighted_null_r{r},95);
                end

            end
        end            
    end   
    % save
    if ~isempty(save_dir)
        save(fullfile(save_dir,[mouse save_suffix]),'-struct','mouse_data','-v7.3')
    end
    disp('     lag data saved')
    
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_exp_dirs: get non-reward directories
function exp_dirs = get_exp_dirs(mouse,data_dir)
    % get a list of experiment directories
    exp_dirs = dir(fullfile(data_dir,mouse));
    is_dirs = [exp_dirs.isdir];
    exp_dirs = {exp_dirs.name}';
    exp_dirs = exp_dirs(is_dirs);
    exp_dirs = exp_dirs(~startsWith(exp_dirs,'.') & ~endsWith(exp_dirs,'r'));
    keep_dirs = ones(size(exp_dirs));
    for d = 1:numel(exp_dirs)
        data = load(fullfile(data_dir,mouse,exp_dirs{d},[mouse '_' exp_dirs{d} '.mat']),'DA','ACh');
        if isempty(data.ACh.Fc_exp_hp_art) || isempty(data.DA.Fc_exp_hp_art)
            keep_dirs(d) = 0;
        end
    end
    exp_dirs = exp_dirs(keep_dirs==1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_session_cross_corrs: get cross corr across sessions
function session_cross_corrs = get_session_cross_corrs(mouse,fib,data_dir,lag)

    % info
    exp_dirs = get_exp_dirs(mouse,data_dir);
    str_rois = fib.ROI_orig(strcmp(fib.mouse,mouse));
    
    % initialize output    
    session_cross_corrs = struct;
    session_cross_corrs.str_rois = str_rois;
    session_cross_corrs.lag = -lag:lag;
    session_cross_corrs.exp_dir = exp_dirs;
    session_cross_corrs.r = nan(numel(str_rois),numel(-lag:lag),numel(exp_dirs));
    session_cross_corrs.n = nan(numel(str_rois),numel(-lag:lag),numel(exp_dirs));
    session_cross_corrs.p = nan(numel(str_rois),numel(-lag:lag),numel(exp_dirs));

    for d = 1:numel(exp_dirs)
        data = load(fullfile(data_dir,mouse,exp_dirs{d},[mouse '_' exp_dirs{d} '.mat']));
        DA = data.DA.Fc_exp_hp_art(:,str_rois);
        ACh = data.ACh.Fc_exp_hp_art(:,str_rois);
        session_cross_corr = get_session_cross_corr(DA,ACh,lag);
        % append (rows = ROIs, cols = lags, slices = sessions)
        session_cross_corrs.r(:,:,d) = session_cross_corr.r;
        session_cross_corrs.n(:,:,d) = session_cross_corr.n;
        session_cross_corrs.p(:,:,d) = session_cross_corr.p;            
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_null_cross_corrs: get null distribution across sessions
function null_cross_corrs = get_null_cross_corrs(mouse,fib,data_dir,lag,varargin)
    %%%  parse optional inputs %%%
    ip = inputParser;
    ip.addParameter('n_it',5000);       % #iterations to run for null distribution
    ip.addParameter('n_timepoints',[]); % #timepoints per iteration (if empty, whole trace)
    ip.parse(varargin{:});
    for j=fields(ip.Results)'
        eval([j{1} '=ip.Results.' j{1} ';']);
    end

    % get info
    exp_dirs = get_exp_dirs(mouse,data_dir);
    str_rois = fib.ROI_orig(strcmp(fib.mouse,mouse));
    
    % initialize output
    null_cross_corrs = struct;
    null_cross_corrs.lag = -lag:lag;            
    null_cross_corrs.exp_dirs = exp_dirs;
    null_cross_corrs.max_abs_r = nan(numel(str_rois),n_it,numel(exp_dirs));
    
    % loop
    for d = 1:numel(exp_dirs)
        % load data
        data = load(fullfile(data_dir,mouse,exp_dirs{d},[mouse '_' exp_dirs{d} '.mat']),'DA','ACh');                        
        DA = data.DA.Fc_exp_hp_art(:,str_rois);
        ACh = data.ACh.Fc_exp_hp_art(:,str_rois);
        disp(['     null: ' exp_dirs{d}])
        null_r = get_session_cross_corr_null(DA,ACh,lag,'n_it',n_it);
        % max statistic permutation: controls family-wise error rate across 
        % lags while respecting the correlation structure between lags
        null_cross_corrs.max_abs_r(:,:,d) = permute(max(abs(null_r),[],2),[1 3 2]);        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_null_session_cross_corr: get null distribution for a single session
function null_r = get_session_cross_corr_null(DA,ACh,lag,varargin)
    %%%  parse optional inputs %%%
    ip = inputParser;    
    ip.addParameter('n_it',5000);       % #iterations to run for null distribution
    ip.addParameter('n_timepoints',[]); % #timepoints per iteration (if empty, whole trace)
    ip.parse(varargin{:});
    for j=fields(ip.Results)'
        eval([j{1} '=ip.Results.' j{1} ';']);
    end
    % preallocate
    null_r = nan(size(DA,2),2*lag+1,n_it);        
    n_timepoints = n_timepoints;
    parfor i = 1:n_it      
        % randomly cut the ACh and switch the positions
        i_cut = randsample((lag+1):(size(DA,1)-lag),1);
        this_DA = DA;
        this_ACh = [ACh((i_cut+1):size(ACh,1),:); ACh(1:i_cut,:)];                        
        if ~isempty(n_timepoints) % if we want a specific number of timepoints
            i_start = randi([lag+1 size(this_DA,1)-n_timepoints-lag]);
            this_DA = this_DA(i_start:(i_start+n_timepoints-1),:);
            this_ACh = this_ACh(i_start:(i_start+n_timepoints-1),:);
        end
        session_cross_corr = get_session_cross_corr(this_DA,this_ACh,lag);
        null_r(:,:,i) = session_cross_corr.r;                
        if rem(i,1000)==0
            disp(['          ' num2str(i)])
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_session_cross_corr: get cross corr for a single session
function session_cross_corr = get_session_cross_corr(DA,ACh,lag)
    
    session_cross_corr = struct;
    session_cross_corr.lag = -lag:lag;
    session_cross_corr.r = nan(size(DA,2),numel(session_cross_corr.lag));
    session_cross_corr.p = nan(size(DA,2),numel(session_cross_corr.lag));
    session_cross_corr.n = nan(size(DA,2),numel(session_cross_corr.lag));
    for r = 1:size(DA,2)
        this_DA = DA(:,r);
        this_ACh = ACh(:,r);
        k = 0;
        for j = -lag:lag
            k = k + 1; % increment for indexing
            lag_ACh = [nan(-min(j,0),1);...
                this_ACh((max(j,0)+1):(end+min(j,0)),1);...
                nan(max(j,0),1)];
            keep_idx = ~isnan(this_DA) & ~isnan(lag_ACh);
            [corr_r,p] = corr(this_DA(keep_idx),lag_ACh(keep_idx));
            session_cross_corr.r(r,k) = corr_r;
            session_cross_corr.p(r,k) = p;
            session_cross_corr.n(r,k) = sum(keep_idx);            
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_session_cross_corr_stats: extract useful across-session stats
% max/min correlations, lags, and significance, based on null        
function output = get_session_cross_corr_stats(session_cross_corrs,alpha_val)

    % some basic characterization and significance testing first
    % note: Fisher-transform r to z via atanh
    mean_r = mean(atanh(session_cross_corrs.r),3); % mean across sessions
    std_r = std(atanh(session_cross_corrs.r),[],3); % std
    null_r = mean(prctile(atanh(session_cross_corrs.null),100-alpha_val,2),3);
    sig_r = abs(mean_r) > null_r; % 
    
    % add some stuff to output: transform back to r for interpretability
    output = struct;
    output.mean = tanh(mean_r);
    output.std = tanh(std_r);
    output.n = size(session_cross_corrs.r,3)*ones(size(mean_r,1),1);
    output.sig = sig_r;
    

    % significant max: require significantly greater than null and is a local max   
    candidate_r = islocalmax(mean_r,2) & sig_r==1 & mean_r>0;
    eligible_local_max = mean_r;
    eligible_local_max(~candidate_r) = nan;
    [max_r,max_lag] = nanmax(eligible_local_max,[],2);
    max_lag = session_cross_corrs.lag(max_lag);
    max_lag(isnan(max_r)) = nan;        
    output.max = [max_r(:) max_lag(:)];

    % significant min: require significantly less than null and is a local min            
    candidate_r = islocalmin(mean_r,2) & sig_r==1 & mean_r<0;
    eligible_local_min = mean_r;
    eligible_local_min(~candidate_r) = nan;
    [min_r,min_lag] = nanmax(eligible_local_min,[],2);
    min_lag = session_cross_corrs.lag(min_lag);
    min_lag(isnan(min_r)) = nan;        
    output.min = [min_r(:) min_lag(:)]; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_required_corr_n: calculate required N (taking into account autocorrelation) 
% to detect desired Pearson correlation strength at specific alpha and power
% see https://sample-size.net/correlation-sample-size/
function n = get_required_corr_n(X,Y,r_target,alpha,power,n_lags)
    
    % first calculate required n required to detect desired Pearson
    % correlation strength at specific alpha and power    
    
    % Fisher z transform of target r
    z_r = atanh(r_target); 
    % Normal quantiles
    z_alpha = norminv(1 - alpha/2);
    z_beta  = norminv(power);    
    % required N
    req_n = ((z_alpha + z_beta) / z_r)^2 + 3; 
    
    % because our signals have autocorrelation, the required n is more than 
    % that; let's figure out what n we need to get that effective n
    n = max([req_n*sum(autocorr(X,n_lags)) req_n*sum(autocorr(Y,n_lags))]);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_cc_lag_sample_stats: extract useful stats from bootstrap data
% max/min correlations, lags, and significance, based on p_value and p_thresh
function output = get_cc_lag_sample_stats(mouse_data,p_thresh,lags)

    output = struct;
    
    % significance
    sig_r = mouse_data.p < p_thresh;

    % significant max: require significance and is a local max   
    candidate_r = islocalmax(mouse_data.r,2) & sig_r==1 & mouse_data.r>0;
    eligible_local_max = mouse_data.r;
    eligible_local_max(~candidate_r) = nan;
    eligible_local_max_lags = repmat(lags,size(mouse_data.r,1),1,size(mouse_data.r,3));
    eligible_local_max_lags(~candidate_r) = nan;   
    for r = 1:size(mouse_data.r,1)
        output.loc_max.r{r} = vec(eligible_local_max(r,:,:));
        output.loc_max.r{r} = output.loc_max.r{r}(~isnan(output.loc_max.r{r}));
        output.loc_max.lags{r} = eligible_local_max_lags(r,:,:);
        output.loc_max.lags{r} = output.loc_max.lags{r}(~isnan(output.loc_max.lags{r}));
    end
    % max
    [max_r,max_lag] = nanmax(eligible_local_max,[],2);
    max_lag = lags(max_lag);
    max_lag(isnan(max_r)) = nan;        
    output.max.r = permute(max_r,[1 3 2]);
    output.max.lag = permute(max_lag,[1 3 2]);

    % significant min: require significance and is a local min            
    candidate_r = islocalmin(mouse_data.r,2) & sig_r==1 & mouse_data.r<0;
    eligible_local_min = mouse_data.r;
    eligible_local_min(~candidate_r) = nan;
    eligible_local_min_lags = repmat(lags,size(mouse_data.r,1),1,size(mouse_data.r,3));
    eligible_local_min_lags(~candidate_r) = nan;   
    for r = 1:size(mouse_data.r,1)        
        output.loc_min.r{r} = vec(eligible_local_min(r,:,:));
        output.loc_min.r{r} = output.loc_min.r{r}(~isnan(output.loc_min.r{r}));
        output.loc_min.lags{r} = eligible_local_min_lags(r,:,:);
        output.loc_min.lags{r} = output.loc_min.lags{r}(~isnan(output.loc_min.lags{r}));        
    end
    
    % min
    [min_r,min_lag] = nanmax(eligible_local_min,[],2);
    min_lag = lags(min_lag);
    min_lag(isnan(min_r)) = nan;        
    output.min.r = permute(min_r,[1 3 2]);
    output.min.lag = permute(min_lag,[1 3 2]);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fit_gmm: fit gaussian mixture model, testing range of #gaussians,
% returning the best on based on goodness of fit
function [gm,gof_vals] = fit_gmm_to_hist(X,varargin)

    %%%  parse optional inputs %%%
    ip = inputParser;
    ip.addParameter('gof','AIC');                   % default AIC; other options 'BIC','NegativeLogLikelihood'
    ip.addParameter('max_n',numel(unique(X(:))));   % the max # of components to run for GMM
    ip.addParameter('choose','min');                % parameter to choose #components; other options: "elbow"
    ip.addParameter('n_tries',3);                   % sometimes the GMM is ill-conditioned: # of retries
    ip.addParameter('MaxIter',1000);                % max #iterations for fitting GMM
    ip.parse(varargin{:});
    for j=fields(ip.Results)'
        eval([j{1} '=ip.Results.' j{1} ';']);
    end
    
    data = X(:);
    
    options = statset('MaxIter',MaxIter);
    gof_vals = nan(max_n,1);
   
    for i = 1:max_n
        n_try = 0;
        while n_try < n_tries
            try
                n_try = n_try + 1;
                gm = fitgmdist(data,i,'Options',options);
                gof_vals(i)= gm.(gof);
                n_try = n_tries;
            catch exception        
            end
        end
    end
    if strcmp(choose,'min')
        [~,n_gauss] = min(gof_vals);
    else
        n_gauss = get_elbow(1:max_n,gof_vals);
    end
    n_try = 0;
    while n_try < n_tries
        try
            n_try = n_try + 1;
            gm = fitgmdist(data,n_gauss,'Options',options);    
            n_try = n_tries;
        catch exception        
        end
    end           
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_elbow
function idx = get_elbow(x,y)

    % get endpoints
    i1 = find(~isnan(y),1,'first');
    i2 = find(~isnan(y),1,'last');
    p1 = [x(i1) y(i1)];
    p2 = [x(i2) y(i2)];

    % distance from points to that line determined by the endpoints
    pts = [vec(x) vec(y)];
    d = point_to_line_distance(pts, p1, p2);

    % get the place of max distance
    if i2==i1 
        idx = i1;
    elseif i2-i1==1
        tmp = [i1 i2];
        [~,idx] = min(y(i1:i2));        
        idx = tmp(idx);
    else
        [~,idx] = max(abs(d));
    end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_gmm_weighted_data: get weighted mean of lagged correlation strength for
% each gaussian in the gaussian mixture model
%
% weighting is a combination of 
%       - component membership via posterior soft weighting 
%       - gaussian weighting according to the mu, sigma of that gaussian

function weighted_data = get_gmm_weighted_mean(data,gm,x)

    % posterior probabilities            
    lag_posteriors = posterior(gm,x(:));
    % gaussian weights
    gauss_weights = cell2mat(arrayfun(@(i) normpdf(x(:),gm.mu(i),gm.Sigma(i)),...
        1:numel(gm.mu),'UniformOutput',false));
    % combined weights
    combined_weights = lag_posteriors.*gauss_weights;
    % weighted data        
    weighted_data = data*combined_weights./sum(combined_weights);    
     
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_main_neg_pos_gmm_idx:
%
% Assumption: at least 2 components: one with ACh-leading (neg), one with 
% DA-leading (pos). Find the most prominent component of each sign

function [comp_idx_neg, comp_idx_pos] = get_main_neg_pos_gmm_idx(gm)

    % first sort by proportion 
    [~,comp_idx] = sort(gm.ComponentProportion,'descend');
    comp_idx_neg = comp_idx(find(gm.mu(comp_idx)<0,1,'first'));
    comp_idx_pos = comp_idx(find(gm.mu(comp_idx)>0,1,'first'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot_lag_hist_and_components:
function f =  plot_lag_hist_and_components(cc_lags_distr,lag,varargin)

    %%%  parse optional inputs %%%
    ip = inputParser;
    ip.addParameter('z_to_r',1);                     % use back-transformed z                       
    ip.addParameter('r_binEdges',[-0.6:0.05:0.6]);   % histogram bin edges
    ip.parse(varargin{:});
    for j=fields(ip.Results)'
        eval([j{1} '=ip.Results.' j{1} ';']);
    end


    % figure
    f = figure('Position',[400 200 550 850]);
    
    % lag distribution
    subplot(3,1,1)
    histogram(cc_lags_distr.lag_hist.gm_mu,'BinEdges',[-18:18],'FaceColor',lines(1),'HandleVisibility','off')
    hold on
    
    for i = 1:numel(cc_lags_distr.lag_gm.all.gm.mu)
        xline(cc_lags_distr.lag_gm.all.gm.mu(i),'-','Color',[.5 .5 .5],'LineWidth',2,'HandleVisibility','off')
    end
    xline(cc_lags_distr.lag_gm.all.gm.mu(cc_lags_distr.lag_gm.all.main_neg_idx),'--b','LineWidth',2)
    xline(cc_lags_distr.lag_gm.all.gm.mu(cc_lags_distr.lag_gm.all.main_pos_idx),'--r','LineWidth',2)
    ylabel('Count')
    xlabel('Lag (frame)')
    set(gca,'YColor',lines(1),'XLim',[-lag lag])
    mus = round([cc_lags_distr.lag_gm.all.gm.mu(cc_lags_distr.lag_gm.all.main_neg_idx),...
        cc_lags_distr.lag_gm.all.gm.mu(cc_lags_distr.lag_gm.all.main_pos_idx)]/lag*1000);
    mus = arrayfun(@(x) num2str(x),mus,'UniformOutput',false);
    mus = cellfun(@(x) [x ' ms'],mus,'UniformOutput',false);
    legend(mus,'Location','best')
    yyaxis right
    plot(-lag:lag,pdf(cc_lags_distr.lag_gm.all.gm,vec(-lag:lag)),'-k','HandleVisibility','off')
    set(gca,'YColor','k')
    ylabel('GMM pdf')
    
    title('Significant Lags')
    
    % negative lag
    subplot(3,1,2)
    hold on
    if z_to_r == 1
        all_mu = tanh(cc_lags_distr.lag_gm.all.weighted_r_z(:,cc_lags_distr.lag_gm.all.main_neg_idx));
        insig_mu = all_mu(cc_lags_distr.lag_gm.all.weighted_r_z_sig(:,cc_lags_distr.lag_gm.all.main_neg_idx)==0);
    else
        all_mu = cc_lags_distr.lag_gm.all.weighted_r(:,cc_lags_distr.lag_gm.all.main_neg_idx);
        insig_mu = all_mu(cc_lags_distr.lag_gm.all.weighted_r_sig(:,cc_lags_distr.lag_gm.all.main_neg_idx)==0);
    end
    histogram(all_mu,'BinEdges',r_binEdges,'FaceColor',[.8 .8 .8],'FaceAlpha',1)
    histogram(insig_mu,'BinEdges',r_binEdges,'FaceColor',[.5 .5 .5],'FaceAlpha',1)
    xlabel('Corr r')
    ylabel('Count')   
    title('ACh \rightarrow DA','interpreter','tex','Color','b')
    
    % positive lag
    subplot(3,1,3)
    hold on
    if z_to_r == 1
        all_mu = cc_lags_distr.lag_gm.all.weighted_r_z(:,cc_lags_distr.lag_gm.all.main_pos_idx);
        insig_mu = all_mu(cc_lags_distr.lag_gm.all.weighted_r_z_sig(:,cc_lags_distr.lag_gm.all.main_pos_idx)==0);
    else
        all_mu = cc_lags_distr.lag_gm.all.weighted_r(:,cc_lags_distr.lag_gm.all.main_pos_idx);
        insig_mu = all_mu(cc_lags_distr.lag_gm.all.weighted_r_sig(:,cc_lags_distr.lag_gm.all.main_pos_idx)==0);
    end    
    histogram(all_mu,'BinEdges',r_binEdges,'FaceColor',[.8 .8 .8],'FaceAlpha',1)
    histogram(insig_mu,'BinEdges',r_binEdges,'FaceColor',[.5 .5 .5],'FaceAlpha',1)
    xlabel('Corr r')
    ylabel('Count')
    title('DA \rightarrow ACh','interpreter','tex','Color','r')

end


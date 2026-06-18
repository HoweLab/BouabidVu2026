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
% 1. cross-correlation across sessions
lag = sr;
cc_results = get_cross_corr_results(mice,fib,data_dir,lag,...
    'percentiles',[2.5 97.5 2.5/(2*lag+1) 97.5/(2*lag+1)],... % including 2 corrected for MC
   'save_by_mouse',1,'save_dir',save_dir1,'n',10000);

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
% results.vals.r = cc_results.r.dominant;
% results.vals.lag = cc_results.lag.dominant;
% results.vals.cluster = vec(kmeans(...
%     [results.vals.r,results.vals.lag],2));
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
% 
% % significant hotspot
% results.sig_moran = local_morans_I_sig_moran(results.smooth,...
%     results.moran,null_moran,'dominant','mask',str.striatum_mask);  
% 
% % save for convenience
% results.str = str;
% save(fullfile(save_dir1,'cross_corr_dominant_results.mat'),'-struct','results')
%     
% %%   
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 4. results figs
% results = load(fullfile(save_dir1,'cross_corr_dominant_results.mat'));
% 
% % scatter plot
% plot_lag_corr_scatter(results.vals.r,results.vals.lag/18*1000,'lat_bins',[-1000:(1000/9):1000]);
% 
% % smooth maps
% outlines = get_mask_projection_outlines(results.sig_moran.map,...
%     results.str,'apply_str_mask',1,'proj_orientations',{'axial','sagittal'});    
% plot_smooth_maps(results.smooth,results.str,'outlines',outlines);
% 
% % in-v-out violin
% plot_violin_in_out(results.vals.r,fib,results.str,results.sig_moran.map)
% 
% % pie chart: mouse composition of hotspot
% map_ind = get_map_ind(fib,results.str);
% in_map = results.sig_moran.map(map_ind) == 1;
% mouse_n = zeros(numel(mice),1);
% hotspot_fib = fib(in_map,:);
% for m = 1:numel(mice)
%     mouse_n(m) =sum(ismember(cellstr(hotspot_fib.mouse),mice{m}));
% end
% figure
% p = pie(mouse_n);
% p_colors = lines(7);
% p_colors(4,:) = p_colors(6,:);
% for i = 1:numel(mice)
%     p((2*i)-1).FaceColor = p_colors(i,:);
% end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_cross_corr_results
function cc_results = get_cross_corr_results(mice,fib,data_dir,lag,varargin)
    
%%%  parse optional inputs %%%
    ip = inputParser;
    ip.addParameter('percentiles',[2.5 97.5]);
    ip.addParameter('n',10000);
    ip.addParameter('save_by_mouse',0);
    ip.addParameter('save_dir',[]);
    ip.parse(varargin{:});
    for j=fields(ip.Results)'
        eval([j{1} '=ip.Results.' j{1} ';']);
    end
    
    % if there's no save dir we don't save anything
    if isempty(save_dir)
        save_by_mouse = 0;
    end
          
    % initialize output
    cc_results = struct;
    cc_results.lag = -lag:lag;
    cc_results.mouse = [];
    cc_results.str_rois = [];
    cc_results.r = [];
    cc_results.sig = [];
    fields_to_fill = {'max_r','max_lag','min_r','min_lag'};
    for f = 1:numel(fields_to_fill)
        cc_results.(fields_to_fill{f}) = [];
    end
    
    for m = 1:numel(mice)        
        mouse = mice{m};
        disp(mouse)        
        if save_by_mouse == 1 && exist(fullfile(save_dir,[mouse '.mat']),'file')
            session_cross_corrs = load(fullfile(save_dir,[mouse '.mat']));
        else
            session_cross_corrs = get_session_cross_corrs(mouse,fib,data_dir,lag);
            null_cross_corrs = get_null_cross_corrs(mouse,fib,data_dir,lag,'percentiles',percentiles,'n',n);
            null_cross_corrs = rmfield(null_cross_corrs,'exp_dirs');
            session_cross_corrs.null = null_cross_corrs;
            if save_by_mouse == 1
                save(fullfile(save_dir,[mouse '.mat']),'-struct','session_cross_corrs','-v7.3')
            end
        end
                        
        % for convenience
        mean_r = mean(session_cross_corrs.r,3);
        r_greater = double(mean_r > mean(session_cross_corrs.null.(['r_'  strrep(num2str(max(percentiles)),'.','p')]),3));
        r_less = -1*double(mean_r < mean(session_cross_corrs.null.(['r_'  strrep(num2str(min(percentiles)),'.','p')]),3));
        
        % add to results                
        cc_results.mouse = [cc_results.mouse; repmat({mouse},numel(session_cross_corrs.str_rois),1)];
        cc_results.str_rois = [cc_results.str_rois; vec(session_cross_corrs.str_rois)];
        cc_results.r = [cc_results.r; mean_r];
        cc_results.sig = [cc_results.sig; r_greater+r_less];        
        
        % some stats
         % max/min correlations, lags, and significance, based on null        
        mean_r = mean(session_cross_corrs.r,3);
       
        % significant max: require significantly greater than null and is a local max
        r_greater = mean_r > mean(session_cross_corrs.null.(['r_'  strrep(num2str(percentiles(p)),'.','p')]),3);
        eligible_local_max = mean_r;
        eligible_local_max(~(islocalmax(mean_r,2) & r_greater)) = nan;
        [session_cross_corrs.max_r,max_lag] = nanmax(eligible_local_max,[],2);
        session_cross_corrs.max_lag = nan(size(max_lag));
        session_cross_corrs.max_lag(~isnan(session_cross_corrs.max_r)) = ...
            cc_results.lag(max_lag(~isnan(session_cross_corrs.max_r)));
        
        % significant min: require significantly less than null and is a local min
        r_less = mean_r < mean(session_cross_corrs.null.(['r_'  strrep(num2str(percentiles(p)),'.','p')]),3);
        eligible_local_min = mean_r;
        eligible_local_min(~(islocalmin(mean_r,2) & r_less)) = nan;
        [session_cross_corrs.min_r,min_lag] = nanmin(eligible_local_min,[],2);        
        session_cross_corrs.min_lag = nan(size(min_lag));
        session_cross_corrs.min_lag(~isnan(session_cross_corrs.min_r)) = ...
            cc_results.lag(min_lag(~isnan(session_cross_corrs.min_r)));
                            
        
        % add to results
        for f = 1:numel(fields_to_fill)
            cc_results.(fields_to_fill{f}) = ...
                [cc_results.(fields_to_fill{f}); vec(session_cross_corrs.(fields_to_fill{f}))];
        end

    end
    
    
    if ~isempty(save_dir)
        save(fullfile(save_dir,'cross_corr_results.mat'),'-struct','cc_results')
    end
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
% get_session_cross_corrs
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
        if ~isempty(data.ACh.Fc_exp_hp_art) && ~isempty(data.DA.Fc_exp_hp_art)
            DA = data.DA.Fc_exp_hp_art(:,str_rois);
            ACh = data.ACh.Fc_exp_hp_art(:,str_rois);
            session_cross_corr = get_session_cross_corr(DA,ACh,lag);
            % append (rows = ROIs, cols = lags, slices = sessions)
            session_cross_corrs.r(:,:,d) = session_cross_corr.r;
            session_cross_corrs.n(:,:,d) = session_cross_corr.n;
            session_cross_corrs.p(:,:,d) = session_cross_corr.p;
            
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_null_cross_corrs
function null_cross_corrs = get_null_cross_corrs(mouse,fib,data_dir,lag,varargin)
    %%%  parse optional inputs %%%
    ip = inputParser;
    ip.addParameter('percentiles',[2.5 97.5]);
    ip.addParameter('n',5000);
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
    for p = 1:numel(percentiles)
        null_cross_corrs.(['r_'  strrep(num2str(percentiles(p)),'.','p')]) = ...
            nan(numel(str_rois),numel(null_cross_corrs.lag),numel(exp_dirs));
    end
    
    % loop
    for d = 1:numel(exp_dirs)
        % load data
        data = load(fullfile(data_dir,mouse,exp_dirs{d},[mouse '_' exp_dirs{d} '.mat']),'DA','ACh');                        
        DA = data.DA.Fc_exp_hp_art(:,str_rois);
        ACh = data.ACh.Fc_exp_hp_art(:,str_rois);
        disp(['     null: ' exp_dirs{d}])
        % preallocate
        session_r = nan(numel(str_rois),numel(null_cross_corrs.lag),n);            
        parfor i = 1:n        
            % randomly cut the ACh and switch the positions
            i_cut = randsample((lag+1):(size(DA,1)-lag),1);
            this_DA = DA;
            this_ACh = [ACh((i_cut+1):size(ACh,1),:); ACh(1:i_cut,:)];                        
            session_cross_corr = get_session_cross_corr(this_DA,this_ACh,lag);                        
            session_r(:,:,i) = session_cross_corr.r;                
            if rem(i,1000)==0
                disp(['          ' num2str(i)])
            end
        end
        for p = 1:numel(percentiles)
            null_cross_corrs.(['r_'  strrep(num2str(percentiles(p)),'.','p')])(:,:,d) = ...
                prctile(session_r,percentiles(p),3);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_session_cross_corr
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
            k = k + 1; % increment for indexint
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

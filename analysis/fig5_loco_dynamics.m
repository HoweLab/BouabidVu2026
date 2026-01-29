%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organization
addpath(fullfile(pwd,'common_functions'))
data_dir = 'D:';
mice = {'UG27','UG28','UG29','UG30','UG31'};
fib = cohort_fib_table(data_dir,mice);
save_dir5 = fullfile(data_dir,'results','5_loco_dynamics');
if ~exist(save_dir5,'dir')
    mkdir(save_dir5)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. concatenate data 
sr = 18; % sampling rate
concat_data = get_concat_data(mice,fib,data_dir,sr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. get initiation & invigoration events
loco_events = get_init_invig_events(concat_data,sr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. get initiation trg avg
init_ta = get_init_trg_avg(concat_data,loco_events,sr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. get invig trg avg (involves fitting GLM to regress out init events)
invig_ta = get_invig_trg_avg(concat_data,loco_events,sr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. save for later
results.init = init_ta;
results.invig = invig_ta;
voxel_size = 0.05;
results.str = get_striatum_vol_mask(...
    [min(fib.fiber_bottom_AP) max(fib.fiber_bottom_AP)],... % AP_range
    [min(fib.fiber_bottom_ML) max(fib.fiber_bottom_ML)],... % ML_range
    [min(fib.fiber_bottom_DV) max(fib.fiber_bottom_DV)],... % DV_range
    voxel_size);
save(fullfile(save_dir5,'loco_trg_avg.mat'),'-struct','results')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6. results figs

% setup
loco_evs = {'init','invig'};
neuromods = {'ACh','DA'};
results = load(fullfile(save_dir5,'loco_trg_avg.mat'));

% pies showing significant dynamics
for i = 1:numel(loco_evs)
    for n = 1:numel(neuromods)
        sig_cols = transpose(struct2array(structfun(@(x) x.([neuromods{n} '_sig'])',...
            results.(loco_evs{i}),'UniformOutput',false)));
        plot_sig_pie(sig_cols)
    end
end

% smooth maps of dominant dynamics
for i = 1:numel(loco_evs)
    for n = 1:numel(neuromods)
        all_ta = struct2array(structfun(@(x) nanmean(x.(neuromods{n}),3),...
            results.(loco_evs{i}),'UniformOutput',false));
        all_ta = all_ta(ismember(-18:36,-9:18),:);
        max_dff = max(all_ta);
        min_dff = min(all_ta);
        dominant_dff = vec(max_dff);
        dominant_dff(abs(min_dff)>abs(max_dff)) = min_dff(abs(min_dff)>abs(max_dff));
        tmp = smooth_activity_map_interp_smoothed(dominant_dff, fib,...
            results.str.info.voxel_size,...
            'AP_range',[min(results.str.info.AP) max(results.str.info.AP)],...
            'ML_range',[min(results.str.info.ML) max(results.str.info.ML)],...
            'DV_range',[min(results.str.info.DV) max(results.str.info.DV)],'gaussian_sigma',2);  
        plot_smooth_maps(tmp.vol_01.smooth,results.str);
    end
end

% smooth maps of triggered averages at various points
timepts_to_plot = -0.25:0.25:1;
sr = 18;
[~,ta_idx] = arrayfun(@(x) min(abs(x - (-18:36))),timepts_to_plot*sr);
for i = 1:numel(loco_evs)
    for n = 1:numel(neuromods)
        all_ta = struct2array(structfun(@(x) nanmean(x.(neuromods{n}),3),...
            results.(loco_evs{i}),'UniformOutput',false));
        all_ta = transpose(all_ta(ta_idx,:));
        tmp = smooth_activity_map_interp_smoothed(all_ta, fib,...
            results.str.info.voxel_size,...
            'AP_range',[min(results.str.info.AP) max(results.str.info.AP)],...
            'ML_range',[min(results.str.info.ML) max(results.str.info.ML)],...
            'DV_range',[min(results.str.info.DV) max(results.str.info.DV)],'gaussian_sigma',2);  
        for j = 1:numel(timepts_to_plot)
            plot_smooth_maps(tmp.(['vol_' sprintf('%02d',j)]).smooth,results.str,...
                'cmap_bounds',[-1 1]*prctile(abs(all_ta(:)),95));
        end
    end
end




%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot_sig_pie
function plot_sig_pie(sig_cols)
% in case i did that thing where i made sig neg column -1
sig_cols = abs(sig_cols); 

% get counts
neither = sum(sig_cols(:,1)==0 & sig_cols(:,2)==0);
pos_only = sum(sig_cols(:,1)==1 & sig_cols(:,2)==0);
neg_only = sum(sig_cols(:,1)==0 & sig_cols(:,2)==1);
both = sum(sig_cols(:,1)==1 & sig_cols(:,2)==1);
pie_vals = [pos_only both neg_only neither];
pie_colors = [1 0 0; .5 0 .5; 0 0 1;.5 .5 .5];

figure
p = pie(pie_vals);
for i = 1:4
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
function concat_data = get_concat_data(mice,fib,data_dir,sr,varargin)

%%%  parse optional inputs %%%
ip = inputParser;
ip.addParameter('nanpad',5*sr)
ip.addParameter('smooth_window',0.3)
ip.addParameter('lp_hz',1.3)
ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end

concat_data = struct;
for m = 1:numel(mice)
    mouse = mice{m};           
    str_rois = fib.ROI_orig(strcmp(fib.mouse,mouse));

    % preallocate
    concat_data.(mouse).date = [];
    concat_data.(mouse).DA = [];
    concat_data.(mouse).ACh = [];
    concat_data.(mouse).vel = [];
    concat_data.(mouse).acc = [];

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
            % append neuromod data
            concat_data.(mouse).DA = [concat_data.(mouse).DA; ...
                nan(nanpad,numel(str_rois)); ...
                data.DA.Fc_exp_hp_art(:,str_rois)];
            concat_data.(mouse).ACh = [concat_data.(mouse).ACh;...
                nan(nanpad,numel(str_rois)); ...
                data.ACh.Fc_exp_hp_art(:,str_rois)];
            % vel & acc
            vel = mean([lowpass(data.behav_ACh.rotaryEncoderVelocity,lp_hz,sr),...
                lowpass(data.behav_DA.rotaryEncoderVelocity,lp_hz,sr)],2);
            vel_smoothed = smooth(vel,sr*smooth_window);
            concat_data.(mouse).vel = [concat_data.(mouse).vel;...
                nan(nanpad,1); vel_smoothed];
            % acc
            acc_smoothed = [smooth(lowpass(diff(vel_smoothed),lp_hz,sr),sr*smooth_window); nan];
            concat_data.(mouse).acc = [concat_data.(mouse).acc;...
                nan(nanpad,1); acc_smoothed];                
            % date
            concat_data.(mouse).date = [concat_data.(mouse).date;...
                repmat({'n/a'},nanpad,1); ...
                repmat(exp_dirs(d),size(data.DA.Fc_exp_hp_art,1),1)];
        end   
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_init_invig_events
function output = get_init_invig_events(concat_data,sr,varargin)

%%%  parse optional inputs %%%
ip = inputParser;
ip.addParameter('still_thresh',0.0208) % calculated from median of mice's velocity while consuming rew
ip.addParameter('still_dur',1) % 1s for stillness
ip.addParameter('loco_dur',2.2) % bout duration calculated from fitting GMM to bout lenths and taking 2nd smallest mean
ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end

mice = fieldnames(concat_data);
for m = 1:numel(mice)
    mouse = mice{m};

    % identify locomotion periods
    loco = double(concat_data.(mouse).vel>still_thresh);
    loco(isnan(concat_data.(mouse).vel)) = nan; % keep appropriate nans
    
    % initiations
    inits = vec(strfind(loco',[zeros(1, still_dur*sr) 1])+still_dur*sr);
    
    % offsets
    all_offsets = strfind(loco',[1 0])+1;
    all_nans = find(isnan(loco));
    
    % find the offsets corresponding to the initiations
    corresponding_offsets = arrayfun(...
        @(x) find(all_offsets>x,1,'first'),inits,'UniformOutput',false);
    if sum(cellfun(@(x) isempty(x),corresponding_offsets))>0 
        corresponding_offsets(cellfun(@(x) isempty(x),corresponding_offsets)) = {nan};
    end
    corresponding_offsets = cell2mat(corresponding_offsets);
    corresponding_offsets(~isnan(corresponding_offsets)) = all_offsets(corresponding_offsets(~isnan(corresponding_offsets)));
    
    % ignore if there are any intervening nans
    corresponding_nans = arrayfun(...
        @(x) find(all_nans>x,1,'first'),inits,'UniformOutput',false);
    if sum(cellfun(@(x) isempty(x),corresponding_nans))>0
        corresponding_nans(cellfun(@(x) isempty(x),corresponding_nans)) = {nan};
    end
    corresponding_nans = cell2mat(corresponding_nans);
    corresponding_nans(~isnan(corresponding_nans)) = all_nans(corresponding_nans(~isnan(corresponding_nans)));
    corresponding_offsets((corresponding_nans-corresponding_offsets)<0) = nan; % if a nan happens before the offset    

    % results
    these_inits = inits(~isnan(corresponding_offsets));
    these_offsets = corresponding_offsets(~isnan(corresponding_offsets));
    output.(mouse).init = these_inits; % initiation
    output.(mouse).term = these_offsets; % termination
 
    % invigorations
    results.(mouse).invig = nan(numel(output.(mouse).init),1);
    for i = 1:numel(output.(mouse).init)
        i1 = output.(mouse).init(i)-1;
        i2 = output.(mouse).term(i)+1;
        this_acc = concat_data.(mouse).acc(i1:i2);
        loc_max = find(islocalmax(this_acc));
        if ~isempty(loc_max)
            results.(mouse).invig(i) = loc_max(1) + i1-1;
        end        
    end 
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get initiation trg avg
function output = get_init_trg_avg(concat_data,loco_events,sr,varargin)

%%%  parse optional inputs %%%
ip = inputParser;
ip.addParameter('loco_dur',2.2) % bout duration calculated from fitting GMM to bout lenths and taking 2nd smallest mean
ip.addParameter('eta_idx',-18:36)
ip.addParameter('idx_of_int',-9:18)
ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end

neuromods = {'ACh','DA'};
mice = fieldnames(concat_data);
for m = 1:numel(mice)
    mouse = mice{m};
    these_inits = loco_events.(mouse).init(...
        (loco_events.(mouse).term-loco_events.(mouse).init)>=(loco_dur*sr));
    for n = 1:numel(neuromods)
        eta = eventTriggeredAverage(concat_data.(mouse).(neuromods{n}), ...
            these_inits, eta_idx(1),eta_idx(end),'nullDistr',1,'bootstrapN',10000);
        eta = eta_sig(eta,'sig_idx',[find(eta_idx==idx_of_int(1)) find(eta_idx==idx_of_int(end))]);
        output.(mouse).(neuromods{n}) = eta.activity;
        output.(mouse).([neuromods{n} '_sig']) = eta.sig_cols;
    end
end
end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get invigoration trg avg
function output = get_invig_trg_avg(concat_data,loco_events,sr,varargin)

%%%  parse optional inputs %%%
ip = inputParser;
ip.addParameter('loco_dur',2.2) % bout duration calculated from fitting GMM to bout lenths and taking 2nd smallest mean
ip.addParameter('eta_idx',-18:36)
ip.addParameter('idx_of_int',-9:18)
ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end

neuromods = {'DA','ACh'};
mice = fieldnames(concat_data);
for m = 1:numel(mice)
    mouse = mice{m};    
    invig_glm = fit_invig_GLM_get_resid(concat_data,loco_events,mouse);
    % now get trg avg on residuals
    for n = 1:numel(neuromods)
        invig_events = find(invig_glm.(neuromods{n}).roi_idx_01.invig_events==1);
        resid_act = struct2array(structfun(@(x) x.invig_resid,invig_glm.(neuromods{n}),...
            'UniformOutput',false));
        invig_sig_pos = structfun(@(x) ~isempty(strfind(transpose(vec(...
            x.mdl.invig.p < 0.05 & x.mdl.invig.mu > 0)),[1 1 1])),invig_glm.(neuromods{n}));
        invig_sig_neg = structfun(@(x) ~isempty(strfind(transpose(vec(...
            x.mdl.invig.p < 0.05 & x.mdl.invig.mu < 0)),[1 1 1])),invig_glm.(neuromods{n}));
        eta = eventTriggeredAverage(resid_act,invig_events,eta_idx(1),eta_idx(end),'nullDistr',0);
        output.(mouse).(neuromods{n}) = eta.activity;
        output.(mouse).([neuromods{n} '_sig']) = [invig_sig_pos invig_sig_neg];
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fit invigoration GLM and get residuals (model out inits, vel, acc)
function output = fit_invig_GLM_get_resid(concat_data,loco_events,mouse,varargin)
    
%%%  parse optional inputs %%%
ip = inputParser;
ip.addParameter('eta_idx',-18:36)
ip.addParameter('idx_of_int',-9:18)
ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end

neuromods = {'ACh','DA'};

data = concat_data.(mouse);

% set up predictors
preds = struct;      
init_events = loco_events.(mouse).init;
invig_events = find(islocalmax(data.acc));
invig_events = invig_events(data.acc(invig_events)>0);

% FIR predictors
for i = 1:numel(idx_of_int)
    % initiation spline        
    preds.(['init' sprintf('%02d',i)]) = zeros(size(data.vel));
    preds.(['init' sprintf('%02d',i)])(init_events + idx_of_int(i)) = 1;

    % invigoration splines         
    preds.(['invig' sprintf('%02d',i)]) = zeros(size(data.vel));
    preds.(['invig' sprintf('%02d',i)])(invig_events + idx_of_int(i)) = 1;        
    preds.(['invig' sprintf('%02d',i)]) = preds.(['invig' sprintf('%02d',i)])(1:numel(data.vel));        
end
preds = orderfields(preds);

% continuous predictors
preds.vel = data.vel;
preds.acc = data.acc;
    
% model specification
pred_names = fieldnames(preds);
continuous_preds = {'vel','acc'};
fir_preds = {'init','invig'};

mdlSpec = ' ~ 1 + vel + acc';
% init
for i = 1:numel(idx_of_int)
    mdlSpec = [mdlSpec ' + init' sprintf('%02d',i)];
end
% invig and interactions
for i = 1:numel(idx_of_int)
    mdlSpec = [mdlSpec ' + invig' sprintf('%02d',i)];
end

% model data around initiations
data_to_model = sort(vec(init_events + eta_idx));

% now fit model
output = struct;
% loop over neuromods
for n = 1:numel(neuromods)
    this_mdlSpec = [neuromods{n} mdlSpec];
    
    % loop over ROIs
    for r = 1:size(data.DA,2)
        
        tbl = struct2table(preds);
        tbl.(neuromods{n}) = data.(neuromods{n})(:,r);
        these_idx = data_to_model(~isnan(tbl.(neuromods{n})(data_to_model)));
        tbl = tbl(these_idx,:);
        
        if ~isempty(tbl)
            mdl = fitglm(tbl,this_mdlSpec);

            tmp = struct;
            tmp.info.idx = these_idx;
            tmp.info.Rsq = mdl.Rsquared.Ordinary;
            tmp.info.Fstat = mdl.devianceTest.FStat;
            tmp.info.df = mdl.DFE;
            tmp.info.AIC = mdl.ModelCriterion.AIC;
            tmp.info.BIC = mdl.ModelCriterion.BIC;
            tmp.info.Coefficients = mdl.Coefficients;
            tmp.pred = mdl.predict;

            % continuous predictors
            for i = 1:numel(continuous_preds)
                tmp.(continuous_preds{i}).mu = ...
                    mdl.Coefficients.Estimate(strcmp(mdl.CoefficientNames,(continuous_preds{i})));
                tmp.(continuous_preds{i}).se = ...
                    mdl.Coefficients.SE(strcmp(mdl.CoefficientNames,(continuous_preds{i})));
                tmp.(continuous_preds{i}).t = ...
                    mdl.Coefficients.tStat(strcmp(mdl.CoefficientNames,(continuous_preds{i})));
                tmp.(continuous_preds{i}).p = ...
                    mdl.Coefficients.pValue(strcmp(mdl.CoefficientNames,(continuous_preds{i})));
            end

            % FIR predictor reconstructions
            for i = 1:numel(fir_preds)
                fir_idx = startsWith(mdl.CoefficientNames,fir_preds{i}) & ~contains(mdl.CoefficientNames,':');
                % now put together
                tmp.(fir_preds{i}).mu = mdl.Coefficients.Estimate(fir_idx,1);
                tmp.(fir_preds{i}).se = mdl.Coefficients.SE(fir_idx,1);
                tmp.(fir_preds{i}).t = mdl.Coefficients.tStat(fir_idx,1);
                tmp.(fir_preds{i}).p = mdl.Coefficients.pValue(fir_idx,1);              
            end
            
            % now get residual using orig table so all ROIs are same size
            % predictors
            pred_mat = table2array(tbl(:,1:end-1)); % ignore the neuromod column
            pred_mat = [ones(size(pred_mat(:,1))) pred_mat]; % add intercept
            % set invig to 0 
            coefs = mdl.Coefficients.Estimate;
            coefs(startsWith(mdl.CoefficientNames,'invig')) = 0;
            % predict 
            pred_sig = pred_mat*coefs;
            % subtract to keep residuals
            resid_sig = tbl.(neuromods{n})-pred_sig;
                   
            % output
            output.(neuromods{n}).(['roi_idx_' sprintf('%02d',r)]).mdl = tmp;
            output.(neuromods{n}).(['roi_idx_' sprintf('%02d',r)]).invig_resid = nan(numel(data_to_model),1);
            output.(neuromods{n}).(['roi_idx_' sprintf('%02d',r)]).invig_resid(these_idx) = vec(resid_sig);
            output.(neuromods{n}).(['roi_idx_' sprintf('%02d',r)]).init_events = tbl.init10;
            output.(neuromods{n}).(['roi_idx_' sprintf('%02d',r)]).invig_events = tbl.invig10;
            
        else
            % fill into our bigger struct
            output.(neuromods{n}).(['roi_idx_' sprintf('%02d',roi_i)]).mdl = [];
            output.(neuromods{n}).(['roi_idx_' sprintf('%02d',roi_i)]).invig_resid = [];
            output.(neuromods{n}).(['roi_idx_' sprintf('%02d',r)]).init_events = tbl.init10;
            output.(neuromods{n}).(['roi_idx_' sprintf('%02d',r)]).invig_events = tbl.invig10;
            clear tmp;
        end
    end    
end
end


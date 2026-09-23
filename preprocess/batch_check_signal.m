%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organization
% addpath(fullfile(pwd,'common_functions'))
addpath('C:\Users\maianhvu\Documents\MATLAB\Scripts\BouabidVu2026\analysis\common_functions');
data_dir = 'G:';
save_dir0 = fullfile(data_dir,'results','0_preprocess');

% mice
mice = struct;
mice.cohort1 = {'UG27','UG28','UG29','UG30','UG31'};
mice.cohort2 = {'AD1','AD2','AD3'};
mice.cohort3 = {'AD4','AD5','AD6'};
mice.cohort4 = {'ADS6','ADS12','ADS13','ADS16','ADS17','ADS19'};
mice.cohort5 = {'609','610','813','816','875','319'};
mice.cohort6 = {'mutAD1','mutAD2','mutAD3'};
mice.cohort7 = {'ADSC_01','ADSC_02','ADSC_03','ADSC_04','ADSC_05'};
mice.cohort8 = {'ADI_02','ADI_03','ADI_04','ADI_05'};
all_mice = struct2cell(structfun(@(x) x(:),mice,'UniformOutput',false));
all_mice = vertcat(all_mice{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. estimate MAD multiplier (k) based on mutant data
% sweep on mutant mice for a duration threshold of 3

% mutant mice and fib
mutant_mice = {'mutAD1','mutAD2','mutAD3'};
fib = cohort_fib_table(data_dir,mutant_mice);

% load if exists
mut_k_results = load_if_exist(fullfile(save_dir0,'mut_k_results.mat'));
% k sweep
if ~isfield(mut_k_results,'mut_k_sweep')
    mut_k_results.k_vals = 1.5:0.5:6;
    mut_k_results.n_val = 3;
    mut_k_results.mut_k_sweep = mutant_k_sweep(data_dir,mutant_mice,fib,'n_val',mut_k_results.n_val,'k_vals',mut_k_results.k_vals);    
    save(fullfile(save_dir0,'mut_k_results.mat'),'-struct','mut_k_results')
end
% estimate a k for each channel, separately for + and -
if ~isfield(mut_k_results,'k_result')
    [~,mut_k_results.k_result] = mutant_k(mut_k_results.mut_k_sweep,mut_k_results.k_vals,'perc_time_thresh',0.05);
    save(fullfile(save_dir0,'mut_k_results.mat'),'-struct','mut_k_results')
end
% all we need is the mad_k
mad_k.green.pos = mut_k_results.k_result.median(1,1);
mad_k.green.neg = mut_k_results.k_result.median(1,2);
mad_k.red.pos = mut_k_results.k_result.median(2,1);
mad_k.red.neg = mut_k_results.k_result.median(2,2);
save(fullfile(save_dir0,'mad_multiplier.mat'),'-struct','mad_k')
clear mut_k_results;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. loop over data and check
sig_check_results = load_if_exist(fullfile(save_dir0,'sig_check_results.mat'));
session_frac_thresh = 0.5; % require a majority of the days to have signal
for m = 1:numel(all_mice)
    mouse = all_mice{m};
    mouse_field = mouse;
    if isstrprop(mouse_field(1),'digit')
        mouse_field = ['m' mouse_field];
    end
    
    % channel names
    if ismember(mouse,mice.cohort1)
        channel_names = {'ACh','DA'};
    elseif ismember(mouse,mice.cohort2)
        channel_names = {'ACh','DA'};
    elseif ismember(mouse,mice.cohort3)
        channel_names = {'AChMut','tdTomato'};
    elseif ismember(mouse,mice.cohort4)
        channel_names = {'ACh','DA'};
    elseif ismember(mouse,mice.cohort5)
        channel_names = {'ACh','DA'};
    elseif ismember(mouse,mice.cohort6)
        channel_names = {'AChMut','DAMut'};
    elseif ismember(mouse,mice.cohort7)
        channel_names = {'AChMut','DA'};
    elseif ismember(mouse,mice.cohort8)
        channel_names = {'ACh','DA'};
    end
    
    if ~isfield(sig_check_results,mouse_field)


        % experiment directories
        exp_dirs = dir(fullfile(data_dir,mouse));
        is_dirs = [exp_dirs.isdir];
        exp_dirs = {exp_dirs.name}';
        exp_dirs = exp_dirs(is_dirs);
        exp_dirs = exp_dirs(~startsWith(exp_dirs,'.') & ~startsWith(exp_dirs,'flagged'));

        % preallocate
        sig_count.(channel_names{1}) = [];
        sig_count.(channel_names{2}) = [];

        % loop
        for d = 1:numel(exp_dirs)
            datapath = fullfile(data_dir,mouse,exp_dirs{d},[mouse '_' exp_dirs{d} '.mat']);
            data = load(datapath);
            if ~isfield(data.(channel_names{1}),'sig') || ~isfield(data.(channel_names{2}),'sig')
                for c = 1:numel(channel_names)
                    channel_name = channel_names{c};            
                    if contains(channel_name,'ACh')
                        pos_mad_k = mad_k.green.pos;
                        neg_mad_k = mad_k.green.neg;
                    else
                        pos_mad_k = mad_k.red.pos;
                        neg_mad_k = mad_k.red.neg;
                    end
                    sig_check = check_data_signal(data.(channel_name),pos_mad_k,neg_mad_k);
                    data.(channel_name).sig = sig_check(:,1);
                end
                save(datapath,'-struct','data','-v7.3')
            end
            % keep running count
            for c = 1:numel(channel_names)
                sig_count.(channel_names{c}) = [sig_count.(channel_names{c}) data.(channel_names{c}).sig];
            end
            disp([mouse ' ' exp_dirs{d}])
        end

        % now see which ones have enough signal        
        
        frac_sess_ch1 = sum(sig_count.(channel_names{1}),2)/size(sig_count.(channel_names{1}),2);
        frac_sess_ch2 = sum(sig_count.(channel_names{2}),2)/size(sig_count.(channel_names{2}),2);
        sig_check_results.(mouse_field) = frac_sess_ch1 > session_frac_thresh & frac_sess_ch2 > session_frac_thresh;
    end    
    save(fullfile(save_dir0,'sig_check_results.mat'),'-struct','sig_check_results');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check_data_signal
% get transients and then use percent time occupied by transients as
% indicator of signal presence
% input is the roi field of the data struct, e.g., data.DA or data.ACh
function sig_check = check_data_signal(data,pos_mad_k,neg_mad_k,varargin)

%%%  parse optional inputs %%%
ip = inputParser;
ip.addParameter('mad_clip_thresh',3);               % clipping bound, in MAD units (separate from detection k)
ip.addParameter('mad_max_iter',10);                 % cap on clipping iterations
ip.addParameter('mad_tol',0.01);                    % relative MAD change to declare convergence
ip.addParameter('tr_n',3);                          % #timepoints in
ip.addParameter('tr_occ_alpha',0.05);               % required percent time occupancy by transients: alpha
ip.addParameter('Fc_field','Fc');                   % which Fc field
ip.addParameter('Fc_artifact_mask','artifact_mask');% artifact mask field; leave blank if n/a
ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end

fc = data.(Fc_field);
if ~isempty(Fc_artifact_mask)
    fc(data.(Fc_artifact_mask)) = nan;
end
transients = get_transients(fc,pos_mad_k,neg_mad_k,...
    'mad_clip_thresh',mad_clip_thresh,'mad_max_iter',mad_max_iter,...
    'mad_tol',mad_tol,'tr_n',tr_n);
tr_perc_time_occ = ...
    vec(sum(transients.transients.sig>0)/size(transients.transients.sig,1) + ... % pos
    sum(transients.transients.sig<0)/size(transients.transients.sig,1));
sig_check = [tr_perc_time_occ>tr_occ_alpha tr_perc_time_occ];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mutant_k
% use the mutant mice k sweep results to get a k multiplier
function [mut_k,k_result] = mutant_k(k_sweep_results,k_vals,varargin)

%%%  parse optional inputs %%%
ip = inputParser;
ip.addParameter('channel_names',{'AChMut','DAMut'});
ip.addParameter('perc_time_thresh',0.05);
ip.addParameter('day_prctile',95);
ip.addParameter('roi_prctile',95);
ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end

% aggregate across days
k_agg_results = aggregate_mutant_k(k_sweep_results);
mice = fieldnames(k_agg_results);
tr_signs = {'pos','neg'};
mut_k = struct;
for m = 1:numel(mice)
    mouse = mice{m};        
    for c = 1:numel(channel_names)        
        channel_name = channel_names{c};
        for s = 1:numel(tr_signs)
            tr_sign = tr_signs{s};
            mut_k.(mouse).(channel_name).(tr_sign) = nan(...                
                size(k_agg_results.(mouse).(channel_name).(tr_sign),2),1);
            for i = 1:size(k_agg_results.(mouse).(channel_name).(tr_sign),2)                
                tmp = find(k_agg_results.(mouse).(channel_name).(tr_sign)(:,i) < perc_time_thresh,1,'first');
                if ~isempty(tmp)
                    mut_k.(mouse).(channel_name).(tr_sign)(i) = k_vals(tmp);
                end                
            end
        end
    end
end

% random sampling w/replacement and bootstrapping 
n_sample = 50;
n_it = 1000;
all_k_results = nan(numel(channel_names),numel(tr_signs),n_it);
for i = 1:n_it
    this_k_result = nan(numel(channel_names),numel(tr_signs),numel(mice)*n_sample);
    for m = 1:numel(mice)
        mouse = mice{m};        
        for c = 1:numel(channel_names)        
            channel_name = channel_names{c};
            for s = 1:numel(tr_signs)
                tr_sign = tr_signs{s};
                tmp = randsample(mut_k.(mouse).(channel_name).(tr_sign),n_sample,1);
                this_k_result(c,s,(m-1)*n_sample+(1:n_sample)) = tmp;
            end
        end
    end
    all_k_results(:,:,i) = prctile(this_k_result,roi_prctile,3);
end
k_result.prctile = prctile(all_k_results,roi_prctile,3);
k_result.median = median(all_k_results,3);
k_result.mean = mean(all_k_results,3);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% aggregate_mutant_k
% aggregate across days via prctile
function mut_k_agg = aggregate_mutant_k(k_sweep_results,varargin)

%%%  parse optional inputs %%%
ip = inputParser;
ip.addParameter('channel_names',{'AChMut','DAMut'});
ip.addParameter('day_prctile',95);
ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end

mut_k_agg = struct;
mice = fieldnames(k_sweep_results);
tr_signs = {'pos','neg'};
for m = 1:numel(mice)
    mouse = mice{m};    
    for c = 1:numel(channel_names)
        channel_name = channel_names{c};
        for s = 1:numel(tr_signs)
            tr_sign = tr_signs{s};
            tmp = prctile(k_sweep_results.(mouse).(channel_name).([tr_sign '_tr_perc_time']),day_prctile,3);
            mut_k_agg.(mouse).(channel_name).(tr_sign) = tmp;
        end
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% k_sweep_mutant: figure out false positive rate with duration threshold
% of n_val (3) samples above/nelow k X MAD amplitude threshold
function mut_k_sweep = mutant_k_sweep(data_dir,mutant_mice,fib,varargin)

%%%  parse optional inputs %%%
ip = inputParser;
ip.addParameter('k_vals',1.5:0.5:5);                % k (multiplier) values to sweep
ip.addParameter('n_val',3);                         % n (duration) value to use
ip.addParameter('sr',18);                           % sr
ip.addParameter('channel_names',{'AChMut','DAMut'});% channels
ip.addParameter('Fc_artifact_mask','artifact_mask');% the field which contains the artifact mask; leave blank if N/A
ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end

mut_k_sweep = struct;

for m = 1:numel(mutant_mice)
    mouse = mutant_mice{m};
    str_rois = fib.ROI_orig(fib.included==1 & ismember(fib.mouse,mouse));
    exp_dirs = get_exp_dirs(mouse,data_dir);
    for c = 1:numel(channel_names)
        channel_name = channel_names{c};
        mut_k_sweep.(mouse).(channel_name).n = nan(numel(exp_dirs),1);
        mut_k_sweep.(mouse).(channel_name).mad = nan(numel(exp_dirs),numel(str_rois));
        mut_k_sweep.(mouse).(channel_name).pos_tr_count = nan(numel(k_vals),numel(str_rois),numel(exp_dirs));
        mut_k_sweep.(mouse).(channel_name).neg_tr_count = nan(numel(k_vals),numel(str_rois),numel(exp_dirs));
        mut_k_sweep.(mouse).(channel_name).pos_tr_hz = nan(numel(k_vals),numel(str_rois),numel(exp_dirs));
        mut_k_sweep.(mouse).(channel_name).neg_tr_hz = nan(numel(k_vals),numel(str_rois),numel(exp_dirs));
        mut_k_sweep.(mouse).(channel_name).pos_tr_perc_time = nan(numel(k_vals),numel(str_rois),numel(exp_dirs));
        mut_k_sweep.(mouse).(channel_name).neg_tr_perc_time = nan(numel(k_vals),numel(str_rois),numel(exp_dirs));
    end 
        
    for d = 1:numel(exp_dirs)
        exp_dir = exp_dirs{d};
        disp([mouse ' ' exp_dir])
        data = load(fullfile(data_dir,mouse,exp_dir,[mouse '_' exp_dir '.mat']));
        for c = 1:numel(channel_names)
            channel_name = channel_names{c};
            % get mad (median absolute deviation)
            if ~isempty(Fc_artifact_mask)
                data.(channel_name).Fc(data.(channel_name).(Fc_artifact_mask)) = nan;
            end
            this_mad = mad(data.(channel_name).Fc(:,str_rois),1,1);
            mut_k_sweep.(mouse).(channel_name).mad(d,:) = this_mad;
            mut_k_sweep.(mouse).(channel_name).n(d) = size(data.(channel_name).Fc,1);
            % now sweep k            
            for k = 1:numel(k_vals)
                k_val = k_vals(k);
                mad_pos = data.(channel_name).Fc(:,str_rois)>k_val*this_mad;
                mad_neg = data.(channel_name).Fc(:,str_rois)<-k_val*this_mad;            
                for r = 1:numel(str_rois)
                    % let's get rid of transients that are too short
                    tmp_pos = mad_pos(:,r)';
                    tmp_neg = mad_neg(:,r)';                    
                    for n = 1:(n_val-1)                        
                        tmp_pos = strrep(tmp_pos,[0 ones(1,n) 0],[0 zeros(1,n) 0]);
                        tmp_neg = strrep(tmp_neg,[0 ones(1,n) 0],[0 zeros(1,n) 0]);
                    end
                    pos_tr_onsets = strfind(tmp_pos,[0 1]);
                    neg_tr_onsets = strfind(tmp_neg,[0 1]);
                    mut_k_sweep.(mouse).(channel_name).pos_tr_count(k,r,d) = numel(pos_tr_onsets);
                    mut_k_sweep.(mouse).(channel_name).neg_tr_count(k,r,d) = numel(neg_tr_onsets);
                    mut_k_sweep.(mouse).(channel_name).pos_tr_hz(k,r,d) = numel(pos_tr_onsets)/(numel(tmp_pos)/sr);
                    mut_k_sweep.(mouse).(channel_name).neg_tr_hz(k,r,d) = numel(neg_tr_onsets)/(numel(tmp_neg)/sr);
                    mut_k_sweep.(mouse).(channel_name).pos_tr_perc_time(k,r,d) = sum(tmp_pos)/sum(~isnan(tmp_pos));
                    mut_k_sweep.(mouse).(channel_name).neg_tr_perc_time(k,r,d) = sum(tmp_neg)/sum(~isnan(tmp_neg));                    
                end                
            end            
        end
    end             
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_exp_dirs: get experiment directories
function exp_dirs = get_exp_dirs(mouse,data_dir,varargin)

    %%%  parse optional inputs %%%
    ip = inputParser;
    ip.addParameter('Fc_field','Fc');                   % which DF/F field to look for
    ip.addParameter('channel_names',{'AChMut','DAMut'});% channel names
    ip.addParameter('datafile_suffix','');              % datafile format [mouse]_[expdir][datafile_suffix].mat; default '', assume format MOUSE_EXPDIR.mat 
    ip.addParameter('excl_unpred_rew',0);               % exclude unpredicted reward recordings (eg UG mice)
    ip.parse(varargin{:});
    for j=fields(ip.Results)'
        eval([j{1} '=ip.Results.' j{1} ';']);
    end
    % get a list of experiment directories
    exp_dirs = dir(fullfile(data_dir,mouse));
    is_dirs = [exp_dirs.isdir];
    exp_dirs = {exp_dirs.name}';
    exp_dirs = exp_dirs(is_dirs);
    exp_dirs = exp_dirs(~startsWith(exp_dirs,'.') & ~startsWith(exp_dirs,'flagged'));
    if excl_unpred_rew == 1
        exp_dirs = exp_dirs(~endsWith(exp_dirs,'r'));
    end
    keep_dirs = ones(size(exp_dirs));
    for d = 1:numel(exp_dirs)      
        data = load(fullfile(data_dir,mouse,exp_dirs{d},[mouse '_' exp_dirs{d} datafile_suffix '.mat']),channel_names{1},channel_names{2});
        if isempty(data.(channel_names{1}).(Fc_field)) || isempty(data.(channel_names{2}).(Fc_field))
            keep_dirs(d) = 0;
        end
    end
    exp_dirs = exp_dirs(keep_dirs==1);
end

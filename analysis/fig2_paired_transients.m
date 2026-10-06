%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organization
addpath(fullfile(pwd,'common_functions'))
data_dir = 'G:';
% the non-pavlovian task mice
coh_mice.cohort1 =  {'UG27','UG28','UG29','UG30','UG31'};
coh_mice.cohort5 = {'609','610','813','816','875'};
mice = struct2cell(structfun(@(x) x(:),coh_mice,'UniformOutput',false));
mice = vertcat(mice{:});

% fib = cohort_fib_table(data_dir,mice,'inclusion_field','in_striatum');
fib = cohort_fib_table(data_dir,mice);
save_dir0 = fullfile(data_dir,'results','0_preprocess');
save_dir1 = fullfile(data_dir,'results','1_cross_corr');

% directory for saving (interim) results
save_dir2 = fullfile(data_dir,'results','2_paired_transients');
if ~exist(save_dir2,'dir')
    mkdir(save_dir2)
end

sr = 18; % sampling rate

% transient detection MAD multiplier
mad_multiplier = load(fullfile(save_dir0,'mad_multiplier.mat'));
mad_multiplier.ACh = mad_multiplier.green;
mad_multiplier.DA = mad_multiplier.red;
mad_multiplier = rmfield(mad_multiplier,{'green','red'});



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. get single transient rate (average across sessions)
tr_stats = load_if_exist(fullfile(save_dir2,'transient_stats.mat'));
if isempty(fieldnames(tr_stats))
    tr_stats = get_transient_stats(fib,data_dir,mad_multiplier,...
        'Fc_field','Fc','channel_names',{'ACh','DA'},'Fc_artifact_mask','artifact_mask',...
        'save_dir',save_dir2);
    save(fullfile(save_dir2,'transient_stats.mat'),'-struct','tr_stats')
end

% 2. get transient info pairing info for all recording locations
paired_tr = load_if_exist(fullfile(save_dir2,'transient_pairing_results.mat'));
if isempty(fieldnames(paired_tr))
    paired_tr = get_all_pairing_stats(fib,data_dir,mad_multiplier,...
        'Fc_field','Fc','channel_names',{'ACh','DA'},'Fc_artifact_mask','artifact_mask',...
        'save_dir',save_dir2,'null_it',1000);
    save(fullfile(save_dir2,'transient_pairing_results.mat'),'-struct','paired_tr')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. maps
% size & striatum mask
map_results = load_if_exist(fullfile(save_dir2,'map_results.mat'));
if ~isfield(map_results,'str')
    tmp = load(fullfile(save_dir1,'map_results'),'str','info');
    map_results.str = tmp.str;
    map_results.info = tmp.info;
    clear tmp;
end
if ~isfield(map_results,'vals')
    map_results.vals = struct;
end

% single stats
trs = fieldnames(tr_stats);
fields_of_int = {'freq','mag'};
for t = 1:numel(trs)
    tr = trs{t};
    for f = 1:numel(fields_of_int)
        field_of_int = fields_of_int{f};
        if ~isfield(map_results.vals,field_of_int)
            % vals
            disp([tr ' ' field_of_int])
            map_results.vals.(tr).(field_of_int) = tr_stats.(tr).(field_of_int);
            % interp
            tmp = get_activity_map_interp(tr_stats.(tr).(field_of_int),fib, map_results.str.info.voxel_size,...
                'AP_range',[min(map_results.str.info.AP) max(map_results.str.info.AP)],...
                'ML_range',[min(map_results.str.info.ML) max(map_results.str.info.ML)],...
                'DV_range',[min(map_results.str.info.DV) max(map_results.str.info.DV)],'incl_plot_info',1);
            map_results.interp.vol.(tr).(field_of_int) = tmp.vol_01.interp;	% interpolated volume (z)
            map_results.interp.n.(tr) = tmp.vol_01.n_mice;   % #mice contrib to each voxel; only need once per
            map_results.interp.F.(tr) = tmp.vol_01.interp_F; % interpolant function; only need once per 
        end
    end
end

% pairing
evs = fieldnames(paired_tr);
fields_of_int = {'lat_mad','lat_mode',...
    'sufficiency','necessity',...
    'sufficiency_corrected','necessity_corrected'};
for e = 1:numel(evs)
    ev = evs{e};
    for f = 1:numel(fields_of_int)
        field_of_int = fields_of_int{f};
        if ~isfield(map_results.vals,field_of_int)
            % vals
            disp([ev ' ' field_of_int])
            map_results.vals.(ev).(field_of_int) = paired_tr.(ev).(field_of_int);
            % interp
            tmp = get_activity_map_interp(paired_tr.(ev).(field_of_int),fib, map_results.str.info.voxel_size,...
                'AP_range',[min(map_results.str.info.AP) max(map_results.str.info.AP)],...
                'ML_range',[min(map_results.str.info.ML) max(map_results.str.info.ML)],...
                'DV_range',[min(map_results.str.info.DV) max(map_results.str.info.DV)],'incl_plot_info',1);
            map_results.interp.vol.(ev).(field_of_int) = tmp.vol_01.interp;	% interpolated volume (z)
            map_results.interp.n.(ev) = tmp.vol_01.n_mice;   % #mice contrib to each voxel; only need once per ev
            map_results.interp.F.(ev) = tmp.vol_01.interp_F; % interpolant function; only need once per ev
        end
    end
end
map_results.info = tmp.info;    
save(fullfile(save_dir2,'map_results.mat'),'-struct','map_results','-v7.3')

% plot smooth maps and fiber maps

        
% % smooth maps and moran calculations
% evs = fieldnames(paired_tr.occurrence);
% for e = 1:numel(evs)
%     ev = evs{e};
%     results = struct;
% 
%     % vals
%     results.vals.occurrence = paired_tr.occurrence.(ev);
%     % add these in for saving convenience
%     results.vals.dir_index = paired_tr.dir_index.(ev);
%     results.vals.mag_corr = paired_tr.mag_corr.(ev)(:,1); 
%     results.vals.mag_corr(results.vals.dir_index<0) = ...
%         paired_tr.mag_corr.(ev)(results.vals.dir_index<0,2);
% 
%     % smooth maps of occurrence
%     tmp = smooth_activity_map_interp_smoothed(results.vals.occurrence, fib,...
%         str.info.voxel_size,'AP_range',[min(str.info.AP) max(str.info.AP)],...
%         'ML_range',[min(str.info.ML) max(str.info.ML)],...
%         'DV_range',[min(str.info.DV) max(str.info.DV)],'gaussian_sigma',2);  
%     results.smooth = tmp.vol_01.smooth;
% 
%     % moran
%     results.moran = local_morans_I(results.smooth,'weight_matrix',ones(21,21,21));
% 
%     % null moran (this can take a while -- easier to run this on the side
%     % and save out results and then come back to it)
%     null_moran = local_morans_I_bootstrap_null(...
%         fib,results.vals.r,voxel_size,...
%         'AP_range',[min(str.info.AP) max(str.info.AP)],...
%         'ML_range',[min(str.info.ML) max(str.info.ML)],...
%         'DV_range',[min(str.info.DV) max(str.info.DV)]);
% 
%     % significant hotspot
%     results.sig_moran = local_morans_I_sig_moran(results.smooth,...
%         results.moran,null_moran,'dominant','mask',str.striatum_mask); 
%     
%     % correlation hotspot comparison
%     results.sig_moran.corr_hotspot_overlap = hotspot_comparison(...
%         results.sig_moran.rand,results.sig_moran.vox,...
%         str.info.DV,str.info.ML,str.info.AP,...
%         corr_hotspot.sig_moran.rand,corr_hotspot.sig_moran.vox,...
%         corr_hotspot.str.info.DV,corr_hotspot.str.info.ML,corr_hotspot.str.info.AP);
% 
%     % save for convenience
%     results.str = str;
%     save(fullfile(save_dir2,[ev '_results.mat']),'-struct','results')
% end
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 3. results figs
% for e = 1:numel(evs)
%     ev = evs{e};
%     
%     % load saved results
%     results = load(fullfile(save_dir2,[ev '_results.mat']));
% 
%     % directionality index histograms
%     plot_directionality_index_histogram(results.vals.dir_index)
%     
%     % smooth maps of occurrence rate, put on contours for hotspots
%     tr_outlines = get_mask_projection_outlines(results.sig_moran.map,...
%         results.str,'apply_str_mask',1,'proj_orientations',{'axial','sagittal'});
%     corr_outlines = get_mask_projection_outlines(corr_hotspot.sig_moran.map,...
%         corr_hotspot.str,'apply_str_mask',1,'proj_orientations',{'axial','sagittal'});
%     outlines.axial = [tr_outlines.axial; corr_outlines.axial];
%     outlines.sagittal = [tr_outlines.sagittal; corr_outlines.sagittal];
%     plot_smooth_maps(results.smooth,results.str,'outlines',outlines);
% 
%     % venn diagrams of hotspot comparisons (title has #voxels)
%     plot_hotspot_comparison_venn(results.sig_moran.corr_hotspot_overlap)
%     
%     % corr hotspot in-v-out violin of mag corr
%     plot_violin_in_out(results.vals.mag_corr,fib,corr_hotspot.str,corr_hotspot.sig_moran.map)
% 
% end



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot_directionality_index_histogram
function plot_directionality_index_histogram(dir_index)
figure
hist_info = histogram(dir_index,'BinEdges',[-1:.1:1],'Normalization','probability');        
hist_vals = hist_info.Values;
hist_edges = hist_info.BinEdges;
bar_centers = .05+(hist_edges(1:end-1));        
cla
hold on
bar(bar_centers(1:numel(hist_vals)/2),hist_vals(1:numel(hist_vals)/2),...
    'BarWidth',1,'FaceColor',[1 0 1])
bar(bar_centers((numel(hist_vals)/2+1):end),hist_vals((numel(hist_vals)/2+1):end),...
    'BarWidth',1,'FaceColor',[0 1 0],'FaceAlpha',.5)
xline(0)
set(gca,'XLim',[-1 1],'XTick',[-1 0 1],'XTickLabel',{'-1 (DA leads)','0','(ACh leads) 1'})
xlabel('Directionality Index')
ylabel('%')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_transient_stats (across mice)
function output = get_transient_stats(fib,data_dir,mad_multiplier,varargin)
ip = inputParser;
ip.addParameter('channel_names',{'ACh','DA'});  % channel names
ip.addParameter('Fc_field','Fc');               % which DF/F field to use
ip.addParameter('Fc_artifact_mask',[]);         % the field which contains the artifact mask; leave blank if N/A
ip.addParameter('sr',18);                       % if supplied, will save interim steps
ip.addParameter('save_dir',[]);                 % if supplied, will save interim steps
ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end
mice = unique(fib.mouse);
signs = {'pos','neg'};

% initialize
output = struct;
for c = 1:numel(channel_names)
    channel_name = channel_names{c};
    for s = 1:numel(signs)
        output.([channel_name '_' signs{s}]).freq = zeros(size(fib,1),1);
        output.([channel_name '_' signs{s}]).mag = zeros(size(fib,1),1);
    end
end

% loop
for m = 1:numel(mice)
    mouse = mice{m};
    mouse_field = get_mouse_field(mouse);
    str_rois = fib.ROI_orig(strcmp(fib.mouse,mouse));
    if ~isempty(save_dir)            
        % transients
        mouse_transients = load_if_exist(fullfile(save_dir,[mouse '_transients.mat']));
        if isempty(fieldnames(mouse_transients))
            disp('   getting transients')
            mouse_transients = concat_transients(mouse,str_rois,data_dir,mad_multiplier, ...
                'Fc_field',Fc_field,'channel_names',channel_names,'Fc_artifact_mask',Fc_artifact_mask,...
                'save_dir',save_dir);
        end
        % get n timepoints and transients for each date
        all_dates = struct2cell(structfun(@(x) unique([...
            x.(channel_names{1}).pos.exp_dir;...
            x.(channel_names{1}).neg.exp_dir;...
            x.(channel_names{2}).pos.exp_dir;...
            x.(channel_names{2}).neg.exp_dir]),...
            mouse_transients,'UniformOutput',false));
        all_dates = unique(vertcat(all_dates{:}));
        all_n = nan(numel(all_dates),1);
        for d = 1:numel(all_dates)
            exp_dir = all_dates{d};
            tp = load(fullfile(data_dir,mouse,exp_dir,[mouse '_' exp_dir '.mat']),[channel_names{1} '_idx']);
            all_n(d) = diff(tp.([channel_names{1} '_idx']))+1;
        end
        roi_fields = fieldnames(mouse_transients);
        for r = 1:numel(roi_fields)
            roi_field = roi_fields{r};
            roi_idx = strcmp(fib.mouse,mouse) & fib.ROI_orig==str2double(roi_field(4:end));
            for c = 1:numel(channel_names)
                channel_name = channel_names{c};
                for s = 1:numel(signs)
                    tr_sign = signs{s};
                    these_tr = mouse_transients.(roi_field).(channel_name).(tr_sign);
                    if ~isempty(these_tr.exp_dir)
                        tr_exp_dirs = unique(these_tr.exp_dir);
                        exp_n = zeros(numel(tr_exp_dirs),2);
                        exp_mag = zeros(numel(tr_exp_dirs),1);
                        for d = 1:numel(tr_exp_dirs)
                            tr_exp_dir = tr_exp_dirs{d};
                            exp_n(d,1) = sum(ismember(these_tr.exp_dir,tr_exp_dir));
                            exp_n(d,2) = all_n(strcmp(all_dates,tr_exp_dir));
                            exp_mag(d) = mean(these_tr.tr_magnitude(ismember(these_tr.exp_dir,tr_exp_dir)));
                        end
                        exp_pr = mean(exp_n(:,1)./exp_n(:,2))*sr;
                        output.([channel_name '_' signs{s}]).freq(roi_idx) = exp_pr;
                        output.([channel_name '_' signs{s}]).mag(roi_idx) = mean(exp_mag);
                    end
                end
            end
        end
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_all_pairing_stats (across mice)
function output = get_all_pairing_stats(fib,data_dir,mad_multiplier,varargin)

ip = inputParser;
ip.addParameter('channel_names',{'ACh','DA'});  % which DF/F field to use
ip.addParameter('Fc_field','Fc');               % which DF/F field to use
ip.addParameter('Fc_artifact_mask',[]);         % the field which contains the artifact mask; leave blank if N/A
ip.addParameter('save_dir',[]);                 % if supplied, will save interim steps
ip.addParameter('null_it',[]);                  % if > 0 will calculate null distribution for pairing stats
ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end

mice = unique(fib.mouse);
signs = {'pos','neg'};

% initialize
output = struct;
for c = 1:numel(channel_names)
    channel_name1 = channel_names{c};
    channel_name2 = channel_names{2/c};
    for s1 = 1:numel(signs)
        sign1 = signs{s1};
        for s2 = 1:numel(signs)
            sign2 = signs{s2};
            ev = [channel_name1 '_' sign1 '_' channel_name2 '_' sign2];
            output.(ev).lat_mad = nan(size(fib,1),1);
            output.(ev).lat_mode = nan(size(fib,1),1);
            output.(ev).sufficiency = nan(size(fib,1),1);
            output.(ev).necessity = nan(size(fib,1),1);
            if null_it > 0
                output.(ev).lat_mad_p = nan(size(fib,1),1);
                output.(ev).sufficiency_pR = nan(size(fib,1),1);
                output.(ev).sufficiency_pL = nan(size(fib,1),1);
                output.(ev).sufficiency_corrected = nan(size(fib,1),1);
                output.(ev).necessity_pR = nan(size(fib,1),1);
                output.(ev).necessity_pL = nan(size(fib,1),1);
                output.(ev).necessity_corrected = nan(size(fib,1),1);
            end
        end
    end
end   
evs = fieldnames(output);

% loop
for m = 1:numel(mice)
    mouse = mice{m};
    mouse_idx = find(ismember(fib.mouse,mouse));
    str_rois = fib.ROI_orig(strcmp(fib.mouse,mouse));
    disp(mouse)
    
    % transient pairing stats
    mouse_pairing_stats = load_if_exist(fullfile(save_dir,[mouse '_transient_pairing_stats.mat']));
    if isempty(fieldnames(mouse_pairing_stats))
        if ~isempty(save_dir)
            
            % transients
            mouse_transients = load_if_exist(fullfile(save_dir,[mouse '_transients.mat']));
            if isempty(fieldnames(mouse_transients))
                disp('   getting transients')
                mouse_transients = concat_transients(mouse,str_rois,data_dir,mad_multiplier, ...
                    'Fc_field',Fc_field,'channel_names',channel_names,'Fc_artifact_mask',Fc_artifact_mask,...
                    'save_dir',save_dir);
            end
            
            % paired transients
            mouse_paired_transients = load_if_exist(fullfile(save_dir,[mouse '_transients_paired.mat']));
            if isempty(fieldnames(mouse_paired_transients))
                disp('   getting transient pairs')
                mouse_paired_transients = get_all_transient_pairs(mouse_transients,...
                    'save_path',fullfile(save_dir,[mouse '_transients_paired.mat']));
            end
            
            % pairing stats
            disp('   getting transient pairing stats')
            mouse_pairing_stats = get_pairing_stats(mouse_transients,mouse_paired_transients,...
                'save_path',fullfile(save_dir,[mouse '_transient_pairing_stats.mat']));            
        else
            disp('   getting transients')
            mouse_transients = concat_transients(mouse,str_rois,data_dir,mad_multiplier, ...
                'Fc_field',Fc_field,'channel_names',channel_names,'Fc_artifact_mask',Fc_artifact_mask);
            disp('   getting transient pairs')
            mouse_paired_transients = get_all_transient_pairs(mouse_transients);
            disp('   getting transient pairing stats')
            mouse_pairing_stats = get_pairing_stats(mouse_transients,mouse_paired_transients);       
            
        end        
    end
    % if we have null distribution data
    if null_it > 0
        % null transient pairing
        mouse_null_paired_transients = load_if_exist(fullfile(save_dir,[mouse '_transients_paired_null.mat']));
        if isempty(fieldnames(mouse_null_paired_transients))
            disp('   getting null transient pairs')
            mouse_null_paired_transients = get_null_transient_pairs(mouse_transients,null_it,...
                'save_path',fullfile(save_dir,[mouse '_transients_paired_null.mat']));
        end
        % null pairing stats
        mouse_null_pairing_stats = load_if_exist(fullfile(save_dir,[mouse '_transient_null_pairing_stats.mat']));
        if isempty(fieldnames(mouse_null_pairing_stats))
            disp('   getting null transient pairing stats')
            mouse_null_pairing_stats = get_null_pairing_stats(mouse_transients,mouse_null_paired_transients,...
                'save_path',fullfile(save_dir,[mouse '_transient_null_pairing_stats.mat']));
        end
    end
    mouse_idx = find(strcmp(fib.mouse,mouse));
    mouse_rois = fib.ROI_orig(mouse_idx);
    for r = 1:numel(mouse_rois)
        roi_idx = mouse_idx(r);
        roi_num = mouse_rois(r);
        roi_field = ['roi' sprintf('%02d',roi_num)];
        for e = 1:numel(evs)
            ev = evs{e};
            if ~isempty(mouse_pairing_stats.(ev).(roi_field).lat_mad)
                output.(ev).lat_mad(roi_idx) = mouse_pairing_stats.(ev).(roi_field).lat_mad;
                output.(ev).lat_mode(roi_idx) = mouse_pairing_stats.(ev).(roi_field).lat_mode;
                output.(ev).sufficiency(roi_idx) = mouse_pairing_stats.(ev).(roi_field).prob(1);
                output.(ev).necessity(roi_idx) = mouse_pairing_stats.(ev).(roi_field).prob(2);
                if null_it > 0
                    output.(ev).lat_mad_p(roi_idx) = sum(...
                        mouse_null_pairing_stats.(ev).(roi_field).lat_mad < ...
                        output.(ev).lat_mad(roi_idx)) / ...
                        numel(mouse_null_pairing_stats.(ev).(roi_field).lat_mad);
                    output.(ev).sufficiency_pR(roi_idx) = sum(...
                        mouse_null_pairing_stats.(ev).(roi_field).prob(:,1) > ...
                        output.(ev).sufficiency(roi_idx)) /...
                        size(mouse_null_pairing_stats.(ev).(roi_field).prob,1);
                    output.(ev).sufficiency_pL(roi_idx) = sum(...
                        mouse_null_pairing_stats.(ev).(roi_field).prob(:,1) < ...
                        output.(ev).sufficiency(roi_idx)) /...
                        size(mouse_null_pairing_stats.(ev).(roi_field).prob,1);
                    output.(ev).sufficiency_corrected(roi_idx) = ...
                        output.(ev).sufficiency(roi_idx) - ...
                        mean(mouse_null_pairing_stats.(ev).(roi_field).prob(:,1));
                    output.(ev).necessity_pR(roi_idx) = sum(...
                        mouse_null_pairing_stats.(ev).(roi_field).prob(:,2) > ...
                        output.(ev).necessity(roi_idx)) /...
                        size(mouse_null_pairing_stats.(ev).(roi_field).prob,1);
                    output.(ev).necessity_pL(roi_idx) = sum(...
                        mouse_null_pairing_stats.(ev).(roi_field).prob(:,2) < ...
                        output.(ev).necessity(roi_idx)) /...
                        size(mouse_null_pairing_stats.(ev).(roi_field).prob,1);
                    output.(ev).necessity_corrected(roi_idx) = ...
                        output.(ev).necessity(roi_idx) - ...
                        mean(mouse_null_pairing_stats.(ev).(roi_field).prob(:,2));
                end
            end
        end
    end     
end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% concat_transients
function output = concat_transients(mouse,mouse_rois,data_dir,mad_multiplier,varargin)

%%%  parse optional inputs %%%
ip = inputParser;
ip.addParameter('channel_names',{'ACh','DA'});      % which DF/F field to use
ip.addParameter('Fc_field','Fc');                   % which DF/F field to use
ip.addParameter('Fc_artifact_mask',[]);             % the field which contains the artifact mask; leave blank if N/A
ip.addParameter('save_dir',[]);                     % if a save directory is supplied, will save
ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end


% some setup
signs = {'pos','neg'}; % peaks and dips

% intialize output
output = struct;
for r = 1:numel(mouse_rois)
    for c = 1:numel(channel_names)
        for s = 1:numel(signs)
            % basic info
            output.(['roi' sprintf('%02d',mouse_rois(r))]).(channel_names{c}).(signs{s}).exp_dir = [];
            output.(['roi' sprintf('%02d',mouse_rois(r))]).(channel_names{c}).(signs{s}).tr_idx = [];
            output.(['roi' sprintf('%02d',mouse_rois(r))]).(channel_names{c}).(signs{s}).tr_onset = [];
            output.(['roi' sprintf('%02d',mouse_rois(r))]).(channel_names{c}).(signs{s}).tr_offset = [];
            output.(['roi' sprintf('%02d',mouse_rois(r))]).(channel_names{c}).(signs{s}).tr_peak = [];
            output.(['roi' sprintf('%02d',mouse_rois(r))]).(channel_names{c}).(signs{s}).tr_magnitude = [];
            output.(['roi' sprintf('%02d',mouse_rois(r))]).(channel_names{c}).(signs{s}).last_rew_del = []; % time since the most recent reward delivery
            output.(['roi' sprintf('%02d',mouse_rois(r))]).(channel_names{c}).(signs{s}).last_rew_con = []; % time since the most recent reward consumption (first lick after deliv)
        end
    end
end


% loop over data
exp_dirs = dir(fullfile(data_dir,mouse));
is_dirs = [exp_dirs.isdir];
exp_dirs = {exp_dirs.name}';
exp_dirs = exp_dirs(is_dirs);
exp_dirs = exp_dirs(startsWith(exp_dirs,'2'));
for d = 1:numel(exp_dirs)
    data = load(fullfile(data_dir,mouse,exp_dirs{d},[mouse '_' exp_dirs{d} '.mat']));
    if ~isempty(data.(channel_names{1}).(Fc_field)) && ~isempty(data.(channel_names{2}).(Fc_field))
        for c = 1:numel(channel_names)
            fc = data.(channel_names{c}).(Fc_field);
            if ~isempty(Fc_artifact_mask)
                fc(data.(channel_names{c}).(Fc_artifact_mask)) = nan;
            end
            tr = get_transients(fc,mad_multiplier.(channel_names{c}).pos,mad_multiplier.(channel_names{c}).neg);
            for s = 1:numel(signs)
                for r = 1:numel(mouse_rois)
                    if data.(channel_names{1}).sig(mouse_rois(r))==1 && data.(channel_names{2}).sig(mouse_rois(r))==1 % signal in both channels

                        % experiment directory
                        output.(['roi' sprintf('%02d',mouse_rois(r))]).(channel_names{c}).(signs{s}).exp_dir = [...
                            output.(['roi' sprintf('%02d',mouse_rois(r))]).(channel_names{c}).(signs{s}).exp_dir;
                            repmat(exp_dirs(d),numel(tr.transients.(signs{s}).peak{mouse_rois(r)}),1)];
                        % give each transient an index
                        output.(['roi' sprintf('%02d',mouse_rois(r))]).(channel_names{c}).(signs{s}).tr_idx = [...
                            output.(['roi' sprintf('%02d',mouse_rois(r))]).(channel_names{c}).(signs{s}).tr_idx;...
                            vec(1:numel(tr.transients.(signs{s}).peak{mouse_rois(r)}))];
                        % idx of onsets
                        output.(['roi' sprintf('%02d',mouse_rois(r))]).(channel_names{c}).(signs{s}).tr_onset = [...
                            output.(['roi' sprintf('%02d',mouse_rois(r))]).(channel_names{c}).(signs{s}).tr_onset;...
                            vec(tr.transients.(signs{s}).onset{mouse_rois(r)})];
                        % idx of offsets
                        output.(['roi' sprintf('%02d',mouse_rois(r))]).(channel_names{c}).(signs{s}).tr_offset = [...
                            output.(['roi' sprintf('%02d',mouse_rois(r))]).(channel_names{c}).(signs{s}).tr_offset;...
                            vec(tr.transients.(signs{s}).offset{mouse_rois(r)})];
                        % idx of peaks (or troughs)
                        output.(['roi' sprintf('%02d',mouse_rois(r))]).(channel_names{c}).(signs{s}).tr_peak = [...
                             output.(['roi' sprintf('%02d',mouse_rois(r))]).(channel_names{c}).(signs{s}).tr_peak;...
                             vec(tr.transients.(signs{s}).peak{mouse_rois(r)})];
                        % DFF magnitude of peaks (or troughs)
                        output.(['roi' sprintf('%02d',mouse_rois(r))]).(channel_names{c}).(signs{s}).tr_magnitude = [...
                             output.(['roi' sprintf('%02d',mouse_rois(r))]).(channel_names{c}).(signs{s}).tr_magnitude;...
                             vec(tr.transients.(signs{s}).magnitude{mouse_rois(r)})];
                        % reward info if relevant (otherwise nan)
                        if endsWith(exp_dirs{d},'r')
                            rew_info = get_reward_info(data.(['behav_' channel_names{c}]));
                            last_rew_del = vec(rew_info.time_since_rew_del(tr.transients.(signs{s}).onset{mouse_rois(r)}));
                            last_rew_con = vec(rew_info.time_since_rew_con(tr.transients.(signs{s}).onset{mouse_rois(r)}));                                  
                        else
                            last_rew_del = nan(numel(tr.transients.(signs{s}).magnitude{mouse_rois(r)}),1);
                            last_rew_con = nan(numel(tr.transients.(signs{s}).magnitude{mouse_rois(r)}),1);
                        end
                        output.(['roi' sprintf('%02d',mouse_rois(r))]).(channel_names{c}).(signs{s}).last_rew_del = [...
                            output.(['roi' sprintf('%02d',mouse_rois(r))]).(channel_names{c}).(signs{s}).last_rew_del;...
                            last_rew_del];
                        output.(['roi' sprintf('%02d',mouse_rois(r))]).(channel_names{c}).(signs{s}).last_rew_con = [...
                            output.(['roi' sprintf('%02d',mouse_rois(r))]).(channel_names{c}).(signs{s}).last_rew_con;...
                            last_rew_con];
                    end
                end
            end
        end
    end
end
if ~isempty(save_dir)
    save(fullfile(save_dir,[mouse '_transients.mat']),'-struct','output')
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_pairing_stats
function output = get_pairing_stats(mouse_transients,mouse_paired_transients,varargin)

%%%  parse optional inputs %%%
ip = inputParser;
ip.addParameter('channel_names',{'ACh','DA'});  % which DF/F field to use
ip.addParameter('rew_ignore',6);                % ignore transients within 6s of rew deliv or consump
ip.addParameter('save_path',[]);                % if supplied, will save pairing stats

ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end

% setup
roi_fields = fieldnames(mouse_transients); 
evs = fieldnames(mouse_paired_transients);

% initialize
output = struct;
for e = 1:numel(evs)
    ev_info = strsplit(evs{e},'_');
    ch1 = ev_info{1};
    sign1 = ev_info{2};
    ch2 = ev_info{3};
    sign2 = ev_info{4};
  
    for r = 1:numel(roi_fields)
        roi_field = roi_fields{r};
        if ~isempty(mouse_paired_transients.(evs{e}).(roi_field))
            % ignore peri-reward period
            filt_ch1 = ~(mouse_transients.(roi_field).(ch1).(sign1).last_rew_del <= rew_ignore | ...
                mouse_transients.(roi_field).(ch1).(sign1).last_rew_con <= rew_ignore);       
            filt_ch2 = ~(mouse_transients.(roi_field).(ch2).(sign2).last_rew_del <= rew_ignore | ...
                mouse_transients.(roi_field).(ch2).(sign2).last_rew_con <= rew_ignore);
        
            % paired transients
            eligible_pairs = mouse_paired_transients.(evs{e}).(roi_field);
            eligible_pairs = eligible_pairs(filt_ch1(eligible_pairs(:,1))==1 & filt_ch2(eligible_pairs(:,2))==1,:);
            output.(evs{e}).(roi_field).tr_idx = eligible_pairs;
            
            % get magnitudes
            output.(evs{e}).(roi_fields{r}).mag = [ ...
                mouse_transients.(roi_field).(ch1).(sign1).tr_magnitude(eligible_pairs(:,1)),...
                mouse_transients.(roi_field).(ch2).(sign2).tr_magnitude(eligible_pairs(:,2))];           
            
            % get latencies
            output.(evs{e}).(roi_fields{r}).lat = ...
                mouse_transients.(roi_field).(ch2).(sign2).tr_peak(eligible_pairs(:,2))-...
                mouse_transients.(roi_field).(ch1).(sign1).tr_peak(eligible_pairs(:,1));
            output.(evs{e}).(roi_field).lat_mad = mad(output.(evs{e}).(roi_field).lat,1);
            output.(evs{e}).(roi_field).lat_mode = mode(output.(evs{e}).(roi_field).lat);
            
            % probabilities
            output.(evs{e}).(roi_fields{r}).prob = [...
                size(eligible_pairs,1) / sum(filt_ch1),...
                size(eligible_pairs,1) / sum(filt_ch2)];                
            
        else
            output.(evs{e}).(roi_field).tr_idx = [];
            output.(evs{e}).(roi_field).mag = [];
            output.(evs{e}).(roi_field).lat = [];
            output.(evs{e}).(roi_field).lat_mad = [];
            output.(evs{e}).(roi_field).prob = [0 0];
        end
    end
    
end
if ~isempty(save_path)
    save(save_path,'-struct','output')
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_null_pairing_stats
function output = get_null_pairing_stats(mouse_transients,mouse_null_paired_transients,varargin)

%%%  parse optional inputs %%%
ip = inputParser;
ip.addParameter('channel_names',{'ACh','DA'});  % which DF/F field to use
ip.addParameter('rew_ignore',6);                % ignore transients within 6s of rew deliv or consump
ip.addParameter('save_path',[]);                % if supplied, will save pairing stats

ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end

% setup
roi_fields = fieldnames(mouse_transients); 
it_fields = fieldnames(mouse_null_paired_transients);
evs = fieldnames(mouse_null_paired_transients.(it_fields{1}));

% initialize
output = struct;

for e = 1:numel(evs)
    ev_info = strsplit(evs{e},'_');
    ch1 = ev_info{1};
    sign1 = ev_info{2};
    ch2 = ev_info{3};
    sign2 = ev_info{4};
    for r = 1:numel(roi_fields)
        roi_field = roi_fields{r};
%         output.(evs{e}).(roi_field).tr_idx = cell(numel(it_fields),1);
%         output.(evs{e}).(roi_field).mag = cell(numel(it_fields),1);
%         output.(evs{e}).(roi_field).lat = cell(numel(it_fields),1);
        output.(evs{e}).(roi_field).lat_mad = nan(numel(it_fields),1);
        output.(evs{e}).(roi_field).lat_mode = nan(numel(it_fields),1);
        output.(evs{e}).(roi_field).prob = nan(numel(it_fields),2);

        for i = 1:numel(it_fields)
            it_field = it_fields{i};

            
            if ~isempty(mouse_null_paired_transients.(it_field).(evs{e}).(roi_field))
                % ignore peri-reward period
                filt_ch1 = ~(mouse_transients.(roi_field).(ch1).(sign1).last_rew_del <= rew_ignore | ...
                    mouse_transients.(roi_field).(ch1).(sign1).last_rew_con <= rew_ignore);       
                filt_ch2 = ~(mouse_transients.(roi_field).(ch2).(sign2).last_rew_del <= rew_ignore | ...
                    mouse_transients.(roi_field).(ch2).(sign2).last_rew_con <= rew_ignore);

                % paired transients
                % note only need to filter on one channel; other was
                % date-scrambled
                eligible_pairs = mouse_null_paired_transients.(it_field).(evs{e}).(roi_field);
                eligible_pairs = eligible_pairs(filt_ch1(eligible_pairs(:,1))==1,:); 
                output.(evs{e}).(roi_field).tr_idx{i} = eligible_pairs;

%                 % get magnitudes
%                 output.(evs{e}).(roi_fields{r}).mag{i} = [ ...
%                     mouse_transients.(roi_field).(ch1).(sign1).tr_magnitude(eligible_pairs(:,1)),...
%                     mouse_transients.(roi_field).(ch2).(sign2).tr_magnitude(eligible_pairs(:,2))];           

                % get latencies
                these_lat = mouse_transients.(roi_field).(ch2).(sign2).tr_peak(eligible_pairs(:,2))-...
                    mouse_transients.(roi_field).(ch1).(sign1).tr_peak(eligible_pairs(:,1));           
%                 output.(evs{e}).(roi_fields{r}).lat{i} = these_lat;
                output.(evs{e}).(roi_field).lat_mad(i) = mad(these_lat,1);
                output.(evs{e}).(roi_field).lat_mode(i) = mode(these_lat);

                % probabilities
                output.(evs{e}).(roi_fields{r}).prob(i,:) = [...
                    size(eligible_pairs,1) / sum(filt_ch1),...
                    size(eligible_pairs,1) / sum(filt_ch2)];                

            end
        end
    end
end
if ~isempty(save_path)
    save(save_path,'-struct','output')
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_a null pairing
function output = get_null_transient_pairs(mouse_transients,n_it,varargin)

%%%  parse optional inputs %%%
ip = inputParser;
ip.addParameter('channel_names',{'ACh','DA'});  % which DF/F field to use
ip.addParameter('rew_ignore',6);                % ignore transients within 6s of rew deliv or consump
ip.addParameter('save_path',[]);                % if supplied, will look for existing file and save results
ip.addParameter('n_update',250);                % how often to display update and save

ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end

% setup
signs = {'pos','neg'};
roi_fields = fieldnames(mouse_transients); 
channel_names = channel_names;
n_update = n_update;

% initialize saved struct if we need to
output = load_if_exist(save_path);
it_to_do = setdiff(arrayfun(@(x) ['it_' sprintf('%05d',x)],vec(1:n_it),'UniformOutput',false),fieldnames(output));
if ~isempty(save_path) && isempty(fieldnames(output)) 
    save(save_path,'-struct','output')
end

% run
if ~isempty(it_to_do)
    disp('   getting null transient pairing')

    for j = 1:ceil(numel(it_to_do)/n_update)
        tmp = cell(n_update,1);
        parfor i = 1:n_update
            this_it = (j-1)*n_update+i;
            disp(this_it)
            
            it_transients = mouse_transients;
            for r = 1:numel(roi_fields)        
                for s = 1:numel(signs)
                    these_exp_dir = it_transients.(roi_fields{r}).(channel_names{2}).(signs{s}).exp_dir;
                    % shuffle the dates on the second channel
                    exp_dates = unique(these_exp_dir);
                    shuffle_dates = exp_dates(randperm(numel(exp_dates)));
                    if ~isempty(these_exp_dir)
                        [~,date_sh_idx] = ismember(these_exp_dir,exp_dates);
                        it_transients.(roi_fields{r}).(channel_names{2}).(signs{s}).exp_dir = shuffle_dates(date_sh_idx);
                    end
                end
            end
            tmp{i} = get_all_transient_pairs(it_transients);
        end
        disp(['      saving ' num2str(j*n_update)])
        % assign and save
        if ~isempty(save_path)
            for i = 1:n_update
                this_it = (j-1)*n_update+i;
                eval(['it_' sprintf('%05d',this_it) ' = tmp{i};'])
                save(save_path,'-append',['it_' sprintf('%05d',this_it)])
            end
        end
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_all_transient_pairs
function output = get_all_transient_pairs(mouse_transients,varargin)

%%%  parse optional inputs %%%
ip = inputParser;
ip.addParameter('channel_names',{'ACh','DA'}); % which DF/F field to use
ip.addParameter('paired_tr_window',18);        % maximum inter-transient spacing
ip.addParameter('save_path',[]);               % if supplied, will save

ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end

% some setup
signs = {'pos','neg'}; % peaks and dips


% loop over pair types
output = struct;
for c = 1:numel(channel_names)
    for s1 = 1:numel(signs)
        for s2 = 1:numel(signs)
            this_pair = [channel_names{c} '_' signs{s1} '_' channel_names{numel(channel_names)/c} '_' signs{s2}];
            output.(this_pair) = get_transient_pairs(mouse_transients,this_pair,paired_tr_window);
        end
    end
end

% save 
if ~isempty(save_path)
    save(save_path,'-struct','output')
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_transient_pairs
function output = get_transient_pairs(mouse_transients,this_pair,paired_tr_window)

% parse event
ev_info = strsplit(this_pair,'_');
tr1 = ev_info(1:2); % first transient event: neuromod + sign
tr2 = ev_info(3:4); % second transient event: neuromod + sign

% loop over ROIs
output = struct;
roi_fields = fieldnames(mouse_transients);
roi_fields = roi_fields(startsWith(roi_fields,'roi'));
for r = 1:numel(roi_fields)
    this_tr1 = mouse_transients.(roi_fields{r}).(tr1{1}).(tr1{2});
    this_tr2 = mouse_transients.(roi_fields{r}).(tr2{1}).(tr2{2});

    % 1. candidate pairs: anchor on second event
    criterion1 = arrayfun(@(x) find(x - this_tr1.tr_onset <= paired_tr_window & x - this_tr1.tr_onset > 0),this_tr2.tr_onset,'UniformOutput',false); % onsets within 1s
    criterion2 = arrayfun(@(x) find(x - this_tr1.tr_peak <= paired_tr_window & x - this_tr1.tr_peak > 0),this_tr2.tr_peak,'UniformOutput',false); % peaks/troughs within 1s
    criterion3 = arrayfun(@(x) find(ismember(this_tr1.exp_dir,x)),this_tr2.exp_dir,'UniformOutput',false); % same recording        
    tmp = cellfun(@intersect,criterion1,criterion2,'UniformOutput',false);
    tmp = cellfun(@intersect,tmp,criterion3,'UniformOutput',false);
    tmp = cellfun(@max,tmp,'UniformOutput',false); % the 1st event closest to the 2nd event
    candidate_pairs1 = [cell2mat(tmp(~cellfun(@isempty,tmp))) find(cellfun(@isempty,tmp)==0)];
    if ~isempty(candidate_pairs1)
        dupl = [0; diff(candidate_pairs1(:,1))==0];
        candidate_pairs1(dupl==1,:) = [];
    end

    % 2. candidate pairs: anchor on first event
    criterion1 = arrayfun(@(x) find(this_tr2.tr_onset - x <= paired_tr_window & this_tr2.tr_onset - x > 0),this_tr1.tr_onset,'UniformOutput',false); % onsets within 1s
    criterion2 = arrayfun(@(x) find(this_tr2.tr_peak - x <= paired_tr_window & this_tr2.tr_peak - x > 0),this_tr1.tr_peak,'UniformOutput',false); % peaks/troughs within 1s
    criterion3 = arrayfun(@(x) find(ismember(this_tr2.exp_dir,x)),this_tr1.exp_dir,'UniformOutput',false); % same recording                
    tmp = cellfun(@intersect,criterion1,criterion2,'UniformOutput',false);
    tmp = cellfun(@intersect,tmp,criterion3,'UniformOutput',false);        
    tmp = cellfun(@min,tmp,'UniformOutput',false); % the 2nd event closest to the 1st event
    candidate_pairs2 = [find(cellfun(@isempty,tmp)==0) cell2mat(tmp(~cellfun(@isempty,tmp)))];
    if ~isempty(candidate_pairs2)
        dupl = [diff(candidate_pairs2(:,2))==0; 0];
        candidate_pairs2(dupl==1,:) = [];
    end

    % 3. candidate pairs: union of 1 & 2, and get rid of duplicates again
    candidate_pairs = union(candidate_pairs1,candidate_pairs2,'rows');
    if ~isempty(candidate_pairs)
        dupl1 = [0; diff(candidate_pairs(:,1))==0];
        candidate_pairs(dupl1==1,:) = [];
        dupl2 = [diff(candidate_pairs(:,2))==0; 0];
        candidate_pairs(dupl2==1,:) = [];    
    end

    % record results (in order of event)
    output.(roi_fields{r}) = candidate_pairs;         
end
end



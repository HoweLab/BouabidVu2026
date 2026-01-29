%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organization
addpath(fullfile(pwd,'common_functions'))
data_dir = 'D:';
mice = {'UG27','UG28','UG29','UG30','UG31'};
fib = cohort_fib_table(data_dir,mice);
save_dir1 = fullfile(data_dir,'results','1_cross_corr');
corr_hotspot = load(fullfile(save_dir1,'cross_corr_dominant_results.mat'),'sig_moran','str');
save_dir2 = fullfile(data_dir,'results','2_paired_transients');
if ~exist(save_dir2,'dir')
    mkdir(save_dir2)
end

sr = 18; % sampling rate



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. get transient pairing info for all recording locations
paired_tr = get_all_pairing_stats(fib,data_dir);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. occurrence hotspot maps

% size & striatum mask
voxel_size = 0.05;
str = get_striatum_vol_mask(...
    [min(fib.fiber_bottom_AP) max(fib.fiber_bottom_AP)],... % AP_range
    [min(fib.fiber_bottom_ML) max(fib.fiber_bottom_ML)],... % ML_range
    [min(fib.fiber_bottom_DV) max(fib.fiber_bottom_DV)],... % DV_range
    voxel_size);

% smooth maps and moran calculations
evs = fieldnames(paired_tr.occurrence);
for e = 1:numel(evs)
    ev = evs{e};
    results = struct;

    % vals
    results.vals.occurrence = paired_tr.occurrence.(ev);
    % add these in for saving convenience
    results.vals.dir_index = paired_tr.dir_index.(ev);
    results.vals.mag_corr = paired_tr.mag_corr.(ev)(:,1); 
    results.vals.mag_corr(results.vals.dir_index<0) = ...
        paired_tr.mag_corr.(ev)(results.vals.dir_index<0,2);

    % smooth maps of occurrence
    tmp = smooth_activity_map_interp_smoothed(results.vals.occurrence, fib,...
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
    save(fullfile(save_dir2,[ev '_results.mat']),'-struct','results')
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. results figs
for e = 1:numel(evs)
    ev = evs{e};
    
    % load saved results
    results = load(fullfile(save_dir2,[ev '_results.mat']));

    % directionality index histograms
    plot_directionality_index_histogram(results.vals.dir_index)
    
    % smooth maps of occurrence rate, put on contours for hotspots
    tr_outlines = get_mask_projection_outlines(results.sig_moran.map,...
        results.str,'apply_str_mask',1,'proj_orientations',{'axial','sagittal'});
    corr_outlines = get_mask_projection_outlines(corr_hotspot.sig_moran.map,...
        corr_hotspot.str,'apply_str_mask',1,'proj_orientations',{'axial','sagittal'});
    outlines.axial = [tr_outlines.axial; corr_outlines.axial];
    outlines.sagittal = [tr_outlines.sagittal; corr_outlines.sagittal];
    plot_smooth_maps(results.smooth,results.str,'outlines',outlines);

    % venn diagrams of hotspot comparisons (title has #voxels)
    plot_hotspot_comparison_venn(results.sig_moran.corr_hotspot_overlap)
    
    % corr hotspot in-v-out violin of mag corr
    plot_violin_in_out(results.vals.mag_corr,fib,corr_hotspot.str,corr_hotspot.sig_moran.map)

end
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
% get_all_pairing_stats (across mice)
function output = get_all_pairing_stats(fib,data_dir)
mice = unique(fib.mouse);
signs = {'pos','neg'};

% initialize
output = struct;
for s1 = 1:numel(signs)
    for s2 = 1:numel(signs)
        ev = ['ACh_' signs{s1} '_DA_' signs{s2}];
        output.dir_index.(ev) = nan(size(fib,1),1);
        output.occurrence.(ev) = nan(size(fib,1),1);
        output.mag_corr.(ev) = nan(size(fib,1),2); % ACh leading, DA leading
    end
end                  
% loop
for m = 1:numel(mice)
    mouse = mice{m};
    mouse_idx = find(ismember(fib.mouse,mouse));
    str_rois = fib.ROI_orig(strcmp(fib.mouse,mouse));
    mouse_transients = concat_transients(mouse,str_rois,data_dir);
    mouse_paired_transients = get_pairing_stats(mouse_transients);
    for s1 = 1:numel(signs)
        for s2 = 1:numel(signs)
            ev = ['ACh_' signs{s1} '_DA_' signs{s2}];
            output.dir_index.(ev)(mouse_idx) = mouse_paired_transients.(ev).dir_index;
            output.occurrence.(ev)(mouse_idx) = mouse_paired_transients.(ev).occ;
            output.mag_corr.(ev)(mouse_idx,:) = mouse_paired_transients.(ev).mag_corr_r;
        end
    end
end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% concat_transients
function output = concat_transients(mouse,mouse_rois,data_dir)

% some setup
signs = {'pos','neg'}; % peaks and dips
neuromods = {'ACh','DA'};

% intialize output
output = struct;
for r = 1:numel(mouse_rois)
    for n = 1:numel(neuromods)
        for s = 1:numel(signs)
            % basic info
            output.(['roi' sprintf('%02d',mouse_rois(r))]).(neuromods{n}).(signs{s}).exp_dir = [];
            output.(['roi' sprintf('%02d',mouse_rois(r))]).(neuromods{n}).(signs{s}).tr_idx = [];
            output.(['roi' sprintf('%02d',mouse_rois(r))]).(neuromods{n}).(signs{s}).tr_onset = [];
            output.(['roi' sprintf('%02d',mouse_rois(r))]).(neuromods{n}).(signs{s}).tr_offset = [];
            output.(['roi' sprintf('%02d',mouse_rois(r))]).(neuromods{n}).(signs{s}).tr_peak = [];
            output.(['roi' sprintf('%02d',mouse_rois(r))]).(neuromods{n}).(signs{s}).tr_magnitude = [];
            output.(['roi' sprintf('%02d',mouse_rois(r))]).(neuromods{n}).(signs{s}).last_rew_del = []; % time since the most recent reward delivery
            output.(['roi' sprintf('%02d',mouse_rois(r))]).(neuromods{n}).(signs{s}).last_rew_con = []; % time since the most recent reward consumption (first lick after deliv)
        end
    end
end


% loop over data
exp_dirs = dir(fullfile(data_dir,mouse));
is_dirs = [exp_dirs.isdir];
exp_dirs = {exp_dirs.name}';
exp_dirs = exp_dirs(is_dirs);
exp_dirs = exp_dirs(~startsWith(exp_dirs,'.'));
for d = 1:numel(exp_dirs)
    data = load(fullfile(data_dir,mouse,exp_dirs{d},[mouse '_' exp_dirs{d} '.mat']));
    if ~isempty(data.ACh.Fc_exp_hp_art) && ~isempty(data.DA.Fc_exp_hp_art)
        for n = 1:numel(neuromods)
            tr = get_transients(data.(neuromods{n}).Fc_exp_hp_art);
            for s = 1:numel(signs)
                for r = 1:numel(mouse_rois)
                    if data.ACh.sig(mouse_rois(r))==1 && data.DA.sig(mouse_rois(r))==1 % signal in both channels

                        % experiment directory
                        output.(['roi' sprintf('%02d',mouse_rois(r))]).(neuromods{n}).(signs{s}).exp_dir = [...
                            output.(['roi' sprintf('%02d',mouse_rois(r))]).(neuromods{n}).(signs{s}).exp_dir;
                            repmat(exp_dirs(d),numel(tr.([signs{s} '_transients']).peak{mouse_rois(r)}),1)];
                        % give each transient an index
                        output.(['roi' sprintf('%02d',mouse_rois(r))]).(neuromods{n}).(signs{s}).tr_idx = [...
                            output.(['roi' sprintf('%02d',mouse_rois(r))]).(neuromods{n}).(signs{s}).tr_idx;...
                            vec(1:numel(tr.([signs{s} '_transients']).peak{mouse_rois(r)}))];
                        % idx of onsets
                        output.(['roi' sprintf('%02d',mouse_rois(r))]).(neuromods{n}).(signs{s}).tr_onset = [...
                            output.(['roi' sprintf('%02d',mouse_rois(r))]).(neuromods{n}).(signs{s}).tr_onset;...
                            vec(tr.([signs{s} '_transients']).onset{mouse_rois(r)})];
                        % idx of offsets
                        output.(['roi' sprintf('%02d',mouse_rois(r))]).(neuromods{n}).(signs{s}).tr_offset = [...
                            output.(['roi' sprintf('%02d',mouse_rois(r))]).(neuromods{n}).(signs{s}).tr_offset;...
                            vec(tr.([signs{s} '_transients']).offset{mouse_rois(r)})];
                        % idx of peaks (or troughs)
                        output.(['roi' sprintf('%02d',mouse_rois(r))]).(neuromods{n}).(signs{s}).tr_peak = [...
                             output.(['roi' sprintf('%02d',mouse_rois(r))]).(neuromods{n}).(signs{s}).tr_peak;...
                             vec(tr.([signs{s} '_transients']).peak{mouse_rois(r)})];
                        % DFF magnitude of peaks (or troughs)
                        output.(['roi' sprintf('%02d',mouse_rois(r))]).(neuromods{n}).(signs{s}).tr_magnitude = [...
                             output.(['roi' sprintf('%02d',mouse_rois(r))]).(neuromods{n}).(signs{s}).tr_magnitude;...
                             vec(tr.([signs{s} '_transients']).magnitude{mouse_rois(r)})];
                        % reward info if relevant (otherwise nan)
                        if endsWith(exp_dirs{d},'r')
                            rew_info = get_reward_info(data.(['behav_' neuromods{n}]));
                            last_rew_del = vec(rew_info.time_since_rew_del(tr.([signs{s} '_transients']).onset{mouse_rois(r)}));
                            last_rew_con = vec(rew_info.time_since_rew_con(tr.([signs{s} '_transients']).onset{mouse_rois(r)}));                                  
                        else
                            last_rew_del = nan(numel(tr.([signs{s} '_transients']).magnitude{mouse_rois(r)}),1);
                            last_rew_con = nan(numel(tr.([signs{s} '_transients']).magnitude{mouse_rois(r)}),1);
                        end
                        output.(['roi' sprintf('%02d',mouse_rois(r))]).(neuromods{n}).(signs{s}).last_rew_del = [...
                            output.(['roi' sprintf('%02d',mouse_rois(r))]).(neuromods{n}).(signs{s}).last_rew_del;...
                            last_rew_del];
                        output.(['roi' sprintf('%02d',mouse_rois(r))]).(neuromods{n}).(signs{s}).last_rew_con = [...
                            output.(['roi' sprintf('%02d',mouse_rois(r))]).(neuromods{n}).(signs{s}).last_rew_con;...
                            last_rew_con];
                    end
                end
            end
        end
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_pairing_stats
function output = get_pairing_stats(mouse_transients)
% setup
rew_ignore = 6; % ignore transients within 6s of rew deliv or consump
roi_fields = fieldnames(mouse_transients); 
mouse_tr_pairs = get_all_transient_pairs(mouse_transients);
neuromods = {'ACh','DA'};
signs = {'pos','neg'};
    
% initialize
output = struct;
for s1 = 1:numel(signs)
    for s2 = 1:numel(signs)
        output.(['ACh_' signs{s1} '_DA_' signs{s2}]).DA_paired = nan(numel(roi_fields),1);
        output.(['ACh_' signs{s1} '_DA_' signs{s2}]).DA_unpaired = nan(numel(roi_fields),1);
        output.(['ACh_' signs{s1} '_DA_' signs{s2}]).ACh_paired = nan(numel(roi_fields),1);
        output.(['ACh_' signs{s1} '_DA_' signs{s2}]).ACh_unpaired = nan(numel(roi_fields),1);
        output.(['ACh_' signs{s1} '_DA_' signs{s2}]).dir_n = nan(numel(roi_fields),2); % ACh-leading, DA-leading  
        output.(['ACh_' signs{s1} '_DA_' signs{s2}]).dir_index = nan(numel(roi_fields),1);
        output.(['ACh_' signs{s1} '_DA_' signs{s2}]).mag_corr_r = nan(numel(roi_fields),2); % ACh-leading, DA-leading        
        output.(['ACh_' signs{s1} '_DA_' signs{s2}]).mag_corr_p = nan(numel(roi_fields),2); % ACh-leading, DA-leading        
        output.(['ACh_' signs{s1} '_DA_' signs{s2}]).peak_lat_mean = nan(numel(roi_fields),2); % latency between peaks
        output.(['ACh_' signs{s1} '_DA_' signs{s2}]).peak_lat_std = nan(numel(roi_fields),2); % latency between peaks
    end
end

evs = fieldnames(output); % paired transient events: note, neuromod order doesn't matter here
% loop over pairs
for e = 1:numel(evs)
    ev_info = strsplit(evs{e},'_');
    ev1 = [ev_info{1} '_' ev_info{2}]; % ACh transient
    ev2 = [ev_info{3} '_' ev_info{4}]; % DA transient
    % loop over ROIs
    for r = 1:numel(roi_fields)
        fib_idx = r;
        filt_ACh = ~(mouse_transients.(roi_fields{r}).(ev_info{1}).(ev_info{2}).last_rew_del <= rew_ignore | ...
            mouse_transients.(roi_fields{r}).(ev_info{1}).(ev_info{2}).last_rew_con <= rew_ignore);       
        filt_DA = ~(mouse_transients.(roi_fields{r}).(ev_info{3}).(ev_info{4}).last_rew_del <= rew_ignore | ...
            mouse_transients.(roi_fields{r}).(ev_info{3}).(ev_info{4}).last_rew_con <= rew_ignore);
        % paired transients
        ACh_first = mouse_tr_pairs.([ev1 '_' ev2]).(roi_fields{r});            
        DA_first = fliplr(mouse_tr_pairs.([ev2 '_' ev1]).(roi_fields{r})); % flip to match order: ACh col 1
        % eligible pairs
        if ~isempty(ACh_first)
            ACh_first = ACh_first(filt_ACh(ACh_first(:,1))==1 & filt_DA(ACh_first(:,2))==1,:);
        end
        if ~isempty(DA_first)
            DA_first = DA_first(filt_ACh(DA_first(:,1))==1 & filt_DA(DA_first(:,2))==1,:);
        end
        % eligible pairs            
        paired_transients = [ACh_first; DA_first];
        if ~isempty(paired_transients)
            output.(evs{e}).ACh_paired(fib_idx) = numel(unique(paired_transients(:,1)));
            output.(evs{e}).DA_paired(fib_idx) = numel(unique(paired_transients(:,2)));
            output.(evs{e}).ACh_unpaired(fib_idx) = sum(filt_ACh)-numel(unique(paired_transients(:,1)));
            output.(evs{e}).DA_unpaired(fib_idx) = sum(filt_DA)-numel(unique(paired_transients(:,2)));           
            % directionality index
            output.(evs{e}).dir_n(fib_idx,:) = [size(ACh_first,1) size(DA_first,1)];
            output.(evs{e}).dir_index(fib_idx) = (size(ACh_first,1)-size(DA_first,1))/(size(ACh_first,1)+size(DA_first,1));
            % magnitude correlation (Pearson): ACh v DA, 
            % ACh-leading pairs
            [corr_r,p] = corr(mouse_transients.(roi_fields{r}).(ev_info{1}).(ev_info{2}).tr_magnitude(ACh_first(:,1)),...
                mouse_transients.(roi_fields{r}).(ev_info{3}).(ev_info{4}).tr_magnitude(ACh_first(:,2)));
            output.(evs{e}).mag_corr_r(fib_idx,1) = corr_r;
            output.(evs{e}).mag_corr_p(fib_idx,1) = p;
            % DA-leading pairs
            [corr_r,p] = corr(mouse_transients.(roi_fields{r}).(ev_info{1}).(ev_info{2}).tr_magnitude(DA_first(:,1)),...
                mouse_transients.(roi_fields{r}).(ev_info{3}).(ev_info{4}).tr_magnitude(DA_first(:,2)));
            output.(evs{e}).mag_corr_r(fib_idx,2) = corr_r;
            output.(evs{e}).mag_corr_p(fib_idx,2) = p;
            % latency: ACh-leading pairs and then DA-leading pairs
            output.(evs{e}).peak_lat_mean(fib_idx,:) = [...
                mean(mouse_transients.(roi_fields{r}).(ev_info{3}).(ev_info{4}).tr_peak(ACh_first(:,2))-...
                mouse_transients.(roi_fields{r}).(ev_info{1}).(ev_info{2}).tr_peak(ACh_first(:,1))),...
                mean(mouse_transients.(roi_fields{r}).(ev_info{1}).(ev_info{2}).tr_peak(DA_first(:,1))-...
                mouse_transients.(roi_fields{r}).(ev_info{3}).(ev_info{4}).tr_peak(DA_first(:,2)))];
            output.(evs{e}).peak_lat_std(fib_idx,:) = [...
                std(mouse_transients.(roi_fields{r}).(ev_info{3}).(ev_info{4}).tr_peak(ACh_first(:,2))-...
                mouse_transients.(roi_fields{r}).(ev_info{1}).(ev_info{2}).tr_peak(ACh_first(:,1))),...
                std(mouse_transients.(roi_fields{r}).(ev_info{1}).(ev_info{2}).tr_peak(DA_first(:,1))-...
                mouse_transients.(roi_fields{r}).(ev_info{3}).(ev_info{4}).tr_peak(DA_first(:,2)))]; 
        end
    end 
    
    % occurrence
    occ_n = [...
        output.(evs{e}).ACh_paired,... % paired ACh
        output.(evs{e}).DA_paired,... % paired DA
        output.(evs{e}).ACh_paired + output.(evs{e}).ACh_unpaired,... % total ACh        
        output.(evs{e}).DA_paired + output.(evs{e}).DA_unpaired]; % total DA
    output.(evs{e}).occ = sum(occ_n(:,1:2),2)./sum(occ_n(:,3:4),2); % occurrence rate: (# transients in pairs) / (total # transients)
    
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_all_transient_pairs
function output = get_all_transient_pairs(mouse_transients)

% some setup
signs = {'pos','neg'}; % peaks and dips
neuromods = {'ACh','DA'};
paired_tr_window = 18; % 1s

% loop over pair types
output = struct;
for n = 1:numel(neuromods)
    for s1 = 1:numel(signs)
        for s2 = 1:numel(signs)
            this_pair = [neuromods{n} '_' signs{s1} '_' neuromods{numel(neuromods)/n} '_' signs{s2}];
            output.(this_pair) = get_transient_pairs(mouse_transients,this_pair,paired_tr_window);
        end
    end
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
    tmp = cellfun(@max,tmp,'UniformOutput',false);
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
    tmp = cellfun(@max,tmp,'UniformOutput',false);
    candidate_pairs2 = [find(cellfun(@isempty,tmp)==0) cell2mat(tmp(~cellfun(@isempty,tmp)))];
    if ~isempty(candidate_pairs2)
        dupl = [diff(candidate_pairs2(:,2))==0; 0];
        candidate_pairs2(dupl==1,:) = [];
    end

    % 3. candidate pairs: union of 1 & 2, and get rid of duplicates again
    candidate_pairs = union(candidate_pairs1,candidate_pairs2,'rows');
    if ~isempty(candidate_pairs)
        dupl1 = [0;diff(candidate_pairs(:,1))==0];
        candidate_pairs(dupl1==1,:) = [];
        dupl2 = [diff(candidate_pairs(:,2))==0; 0];
        candidate_pairs(dupl2==1,:) = [];    
    end

    % record results (in order of event)
    output.(roi_fields{r}) = candidate_pairs;         
end
end



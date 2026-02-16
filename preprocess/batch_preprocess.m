% preprocess cohorts
% roi1 and behav1 correspond to ttlIn1, green camera (470nm, ACh)
% roi2 and behav2 correspond to ttlIn2, red camera (470nm, ACh)

data_dir = 'D:';

% cohort1
mice = {'UG27','UG28','UG29','UG30','UG31'};
for m = 1:numel(mice)
    mouse = mice{m};
    exp_dirs = dir(fullfile(data_dir,mouse));
    is_dirs = [exp_dirs.isdir];
    exp_dirs = {exp_dirs.name}';
    exp_dirs = exp_dirs(is_dirs);
    exp_dirs = exp_dirs(~startsWith(exp_dirs,'.'));
    for d = 1:numel(exp_dirs)
        disp([mouse ' ' exp_dirs{d}])
        savepath = fullfile(data_dir,mouse,exp_dirs{d},[mouse '_' exp_dirs{d} '.mat']);        
        if ~exist(savepath,'file')
            files = dir(fullfile(data_dir,mouse,exp_dirs{d}));
            files = {files.name}';
            % roi files
            roi_files = files(contains(files,'ROIs') & endsWith(files,'.mat'));
            path_roi1 = fullfile(data_dir,mouse,exp_dirs{d},roi_files{~startsWith(roi_files,'R')});
            path_roi2 = fullfile(data_dir,mouse,exp_dirs{d},roi_files{startsWith(roi_files,'R')});
            % beh files
            beh_files = files(contains(files,'ttlIn') & contains(files,'movie') & endsWith(files,'.mat'));
            path_behav1 = fullfile(data_dir,mouse,exp_dirs{d},beh_files{contains(beh_files,'ttlIn1')});
            path_behav2 = fullfile(data_dir,mouse,exp_dirs{d},beh_files{contains(beh_files,'ttlIn2')});
            % artifacts
            if ismember('artifacts.mat',files)
                art_info = load(fullfile(data_dir,mouse,exp_dirs{d},'artifacts.mat'));
                rm_neg_artifact = [art_info.ACh_neg art_info.DA_neg];
                rm_pos_artifact = [art_info.ACh_pos art_info.DA_pos];
            else
                rm_neg_artifact = [0 0];
                rm_pos_artifact = [0 0];
            end
            % preprocess and save
            data = preprocess_data(path_roi1,path_roi2,path_behav1,path_behav2,...
                'rm_neg_artifact',rm_neg_artifact,'rm_pos_artifact',rm_pos_artifact);  
            output = struct;
            output.ACh = data.roi1;           
            output.DA = data.roi2;
            % in case we didn't need artifact correction, just put this one
            % there so we can use the same field for analyses going forward
            if ~isfield(output.ACh,'Fc_exp_hp_art')
                output.ACh.Fc_exp_hp_art = output.ACh.Fc_exp_hp;
            end
            if ~isfield(output.DA,'Fc_exp_hp_art')
                output.DA.Fc_exp_hp_art = output.DA.Fc_exp_hp;
            end
            output.behav_ACh = data.behav1;
            output.behav_DA = data.behav2;
            output.ACh_idx = data.roi1idx;
            output.DA_idx = data.roi2idx;
            % whether we have signal that day
            output.ACh.sig = 0;
            if ~isempty(output.ACh.Fc_exp_hp_art)
                output.ACh.sig = get_signal_tr(output.ACh.Fc_exp_hp_art);            
            end
            output.DA.sig = 0;
            if ~isempty(output.DA.Fc_exp_hp_art)
                output.DA.sig = get_signal_tr(output.DA.Fc_exp_hp_art);            
            end
            save(savepath,'-struct','output')
        end
    end
end
   
% cohort 2
mice = {'AD1','AD2','AD3','AD4','AD5','AD6'};
for m = 1:numel(mice)
    mouse = mice{m};
    exp_dirs = dir(fullfile(data_dir,mouse));
    is_dirs = [exp_dirs.isdir];
    exp_dirs = {exp_dirs.name}';
    exp_dirs = exp_dirs(is_dirs);
    exp_dirs = exp_dirs(~startsWith(exp_dirs,'.'));
    for d = 1:numel(exp_dirs)
        disp([mouse ' ' exp_dirs{d}])
        savepath = fullfile(data_dir,mouse,exp_dirs{d},[mouse '_' exp_dirs{d} '.mat']);
        if ~exist(savepath,'file')
            files = dir(fullfile(data_dir,mouse,exp_dirs{d}));
            files = {files.name}';
            % roi files
            roi_files = files(contains(files,'ROIs') & endsWith(files,'.mat'));
            path_roi1 = fullfile(data_dir,mouse,exp_dirs{d},roi_files{contains(roi_files,'470')});
            path_roi2 = fullfile(data_dir,mouse,exp_dirs{d},roi_files{contains(roi_files,'570')});
            % beh files
            beh_files = files(contains(files,'ttlIn') & contains(files,'movie') & endsWith(files,'.mat'));
            path_behav1 = fullfile(data_dir,mouse,exp_dirs{d},beh_files{contains(beh_files,'ttlIn1')});
            path_behav2 = fullfile(data_dir,mouse,exp_dirs{d},beh_files{contains(beh_files,'ttlIn2')});
            task_name = strsplit(beh_files{contains(beh_files,'ttlIn1')},'_');
            task_name = task_name{2};
            % preprocess and save
            data = preprocess_data(path_roi1,path_roi2,path_behav1,path_behav2);  
            output = struct;
            output.ACh = data.roi1;           
            output.DA = data.roi2;
            output.behav_ACh = data.behav1;
            output.behav_DA = data.behav2;
            output.ACh_idx = data.roi1idx;
            output.DA_idx = data.roi2idx;  
            % whether we have signal that day
            output.ACh.sig = 0;
            if ~isempty(output.ACh.Fc_exp_hp)
                output.ACh.sig = get_signal_tr(output.ACh.Fc_exp_hp);            
            end
            output.DA.sig = 0;
            if ~isempty(output.DA.Fc_exp_hp)
                output.DA.sig = get_signal_tr(output.DA.Fc_exp_hp);            
            end
            output.task = task_name;
            save(savepath,'-struct','output')
        end
    end
end


% cohort 3: we don't want the high-pass filter on this one so just ignore
% that output
mice = {'ADS6','ADS12','ADS13','ADS16','ADS17','ADS19'};
for m = 1:numel(mice)
    mouse = mice{m};
    exp_dirs = dir(fullfile(data_dir,mouse));
    is_dirs = [exp_dirs.isdir];
    exp_dirs = {exp_dirs.name}';
    exp_dirs = exp_dirs(is_dirs);
    exp_dirs = exp_dirs(~startsWith(exp_dirs,'.'));
    for d = 1:numel(exp_dirs)
        disp([mouse ' ' exp_dirs{d}])
        savepath = fullfile(data_dir,mouse,exp_dirs{d},[mouse '_' exp_dirs{d} '.mat']);        
        if ~exist(savepath,'file')
            files = dir(fullfile(data_dir,mouse,exp_dirs{d}));
            files = {files.name}';
            % roi files
            roi_files = files(contains(files,'ROIs') & endsWith(files,'.mat'));
            path_roi1 = fullfile(data_dir,mouse,exp_dirs{d},roi_files{startsWith(roi_files,'G')});
            path_roi2 = fullfile(data_dir,mouse,exp_dirs{d},roi_files{startsWith(roi_files,'R')});
            % beh files
            beh_files = files(contains(files,'ttlIn') & contains(files,'movie') & endsWith(files,'.mat'));
            path_behav1 = fullfile(data_dir,mouse,exp_dirs{d},beh_files{contains(beh_files,'ttlIn1')});
            path_behav2 = fullfile(data_dir,mouse,exp_dirs{d},beh_files{contains(beh_files,'ttlIn2')});
            task_name = strsplit(beh_files{contains(beh_files,'ttlIn1')},'_');
            task_name = task_name{2};
            % preprocess and save
            data = preprocess_data(path_roi1,path_roi2,path_behav1,path_behav2);  
            output = struct;
            output.ACh = data.roi1;
            output.DA = data.roi2;
            output.behav_ACh = data.behav1;
            output.behav_DA = data.behav2;
            output.ACh_idx = data.roi1idx;
            output.DA_idx = data.roi2idx;   
            % whether we have signal that day
            output.ACh.sig = 0;
            if ~isempty(output.ACh.Fc_exp)
                output.ACh.sig = get_signal_tr(output.ACh.Fc_exp);            
            end
            output.DA.sig = 0;
            if ~isempty(output.DA.Fc_exp)
                output.DA.sig = get_signal_tr(output.DA.Fc_exp);            
            end
            save(savepath,'-struct','output')
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% has signal?
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sig = get_signal_tr(Fc)
    transients = get_transients(Fc);   
    n_peaks = cellfun(@(x) numel(x),transients.pos_transients.peak);
    n_troughs = cellfun(@(x) numel(x),transients.neg_transients.peak);
    % we'll assume each transient is minimum duration (3 timepts), though
    % many are longer
    sig_peaks = 3*100*n_peaks./transpose(sum(~isnan(transients.roi)))>2.5;
    sig_troughs = 3*100*n_troughs./transpose(sum(~isnan(transients.roi)))>2.5;
    % we have signal if the rate of detected peaks or detected troughs 
    % exceeds chance
    sig = sig_peaks | sig_troughs;
end
       
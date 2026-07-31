% preprocess cohorts
% roi1 and behav1 correspond to ttlIn1, green camera (470nm, ACh)
% roi2 and behav2 correspond to ttlIn2, red camera (470nm, ACh)

data_dir = 'G:';
mice.cohort1 =  {'UG27','UG28','UG29','UG30','UG31'};
mice.cohort2 = {'AD1','AD2','AD3'};
mice.cohort3 = {'AD4','AD5','AD6'};
mice.cohort4 = {'ADS6','ADS12','ADS13','ADS16','ADS17','ADS19'};
all_mice = struct2cell(structfun(@(x) x(:),mice,'UniformOutput',false));
all_mice = vertcat(all_mice{:});

for m = 1:numel(all_mice)
    mouse = all_mice{m};
    exp_dirs = dir(fullfile(data_dir,mouse));
    is_dirs = [exp_dirs.isdir];
    exp_dirs = {exp_dirs.name}';
    exp_dirs = exp_dirs(is_dirs);
    exp_dirs = exp_dirs(~startsWith(exp_dirs,'.'));
    for d = 1:numel(exp_dirs)
        disp([mouse ' ' exp_dirs{d}])
        savepath = fullfile(data_dir,mouse,exp_dirs{d},[mouse '_' exp_dirs{d} '.mat']);        
%         if ~exist(savepath,'file')
            files = dir(fullfile(data_dir,mouse,exp_dirs{d},'intermed'));
            files = {files.name}';
            % beh files
            beh_files = files(contains(files,'ttlIn') & contains(files,'movie') & endsWith(files,'.mat'));
            path_behav1 = fullfile(data_dir,mouse,exp_dirs{d},'intermed',beh_files{contains(beh_files,'ttlIn1')});
            path_behav2 = fullfile(data_dir,mouse,exp_dirs{d},'intermed',beh_files{contains(beh_files,'ttlIn2')}); 
            % roi files & other cohort specific settings
            roi_files = files(contains(files,'ROIs') & endsWith(files,'.mat'));
            if ismember(mouse,mice.cohort1)
                path_roi1 = fullfile(data_dir,mouse,exp_dirs{d},'intermed',roi_files{~startsWith(roi_files,'R')});
                path_roi2 = fullfile(data_dir,mouse,exp_dirs{d},'intermed',roi_files{startsWith(roi_files,'R')});
                channel_names = {'ACh','DA'};
                light_blink_artifact = [0 1];
            elseif ismember(mouse,mice.cohort2)
                path_roi1 = fullfile(data_dir,mouse,exp_dirs{d},'intermed',roi_files{contains(roi_files,'470')});
                path_roi2 = fullfile(data_dir,mouse,exp_dirs{d},'intermed',roi_files{contains(roi_files,'570')});
                channel_names = {'ACh','DA'};
                light_blink_artifact = [0 0];
            elseif ismember(mouse,mice.cohort3)
                path_roi1 = fullfile(data_dir,mouse,exp_dirs{d},'intermed',roi_files{contains(roi_files,'470')});
                path_roi2 = fullfile(data_dir,mouse,exp_dirs{d},'intermed',roi_files{contains(roi_files,'570')});
                channel_names = {'AChMut','tdTomato'};
                light_blink_artifact = [0 0];                
            elseif ismember(mouse,mice.cohort4)
                path_roi1 = fullfile(data_dir,mouse,exp_dirs{d},'intermed',roi_files{startsWith(roi_files,'G')});
                path_roi2 = fullfile(data_dir,mouse,exp_dirs{d},'intermed',roi_files{startsWith(roi_files,'R')});
                channel_names = {'ACh','DA'};
                light_blink_artifact = [0 0];
            end
            
            % preprocess and save
            output = preprocess_data(path_roi1,path_roi2,path_behav1,path_behav2,...
                'channel_names',channel_names, 'light_blink_artifact',light_blink_artifact); 
            
            % task name for mice performing the pavlovian task
            if ismember(mouse,mice.cohort2) || ismember(mouse,mice.cohort3)
                task_name = strsplit(beh_files{contains(beh_files,'ttlIn1')},'_');
                task_name = task_name{2};
                output.task = task_name;
            end
            
            % save
            save(savepath,'-struct','output')
%         end
    end
end
   



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organization
addpath(fullfile(pwd,'common_functions'))
data_dir = 'D:';
mice = {'UG27','UG28','UG29','UG30','UG31'};
fib = cohort_fib_table(data_dir,mice);
save_dir1 = fullfile(data_dir,'results','1_cross_corr');
corr_hotspot = load(fullfile(save_dir1,'cross_corr_dominant_results.mat'),'sig_moran','str');
save_dir3 = fullfile(data_dir,'results','3_unpred_rew');
if ~exist(save_dir3,'dir')
    mkdir(save_dir3)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 1. concatenate data 
nanpad = sr*5; % put a bunch of nans between days so the cross-correlation sliding doesn't combine data across days
concat_data = get_concat_data(mice,fib,data_dir,nanpad);

% Mai-Anh Vu, updated 2025/08/29
function eta = eta_combine(eta_cell,varargin)
%
% eta_cell should be a nx2 cell array, with the first column being the eta
% file paths or structs, and the second column being the data file paths or 
% data matrix itself (this function will replace roi struct with roi.Fc
%
% This assumes you're combining event triggered average structs with
% similar parameters (windowIdx, etc) and that the only difference is which
% recording you're pulling from.
%
%

%%%  parse optional inputs %%%
ip = inputParser;
ip.addParameter('bootstrap_n',10000);
ip.addParameter('n_rounds',10);
ip.addParameter('roi_nums',[]);
ip.addParameter('hp',[]); % highpass (empty [] for nothing, otherwise [freq sampling_rate])
ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end
if n_rounds < 2
    n_rounds = 2;
end

%%% load where necessary, also make note of bootstrapping
if bootstrap_n > 0
    yes_bootstrap = 1;
    for i = 1:size(eta_cell,1)
        for j = 1:size(eta_cell,2)
            if ischar(eta_cell{i,j})
                if ~endsWith(eta_cell{i,j},'.mat')
                    eta_cell{i,j} = [eta_cell{i,j} '.mat'];
                end
                if j == 1 % only load the eta structs
                    eta_cell{i,j} = load(eta_cell{i,1});
                end
            end    
        end

        % if the second one is an roi struct, replace it with roi.Fc
        if isstruct(eta_cell{i,2})
            if isfield(eta_cell{i,2},'Fc')
                eta_cell{i,2} = eta_cell{i,2}.Fc;
            end
        end
        if ~isfield(eta_cell{i,1},'bootstrap')
            yes_bootstrap = 0;
        end
    end
else
    yes_bootstrap = 0;
end

%%% let's combine some things
eta = struct;
eta.activity = [];
eta.windowIdx = [];
for i = 1:size(eta_cell,1)
    eta.activity = cat(3,eta.activity,eta_cell{i,1}.activity);
    eta.events{i} = eta_cell{i,1}.events;
    eta.windowIdx = cat(2,eta.windowIdx,eta_cell{i,1}.windowIdx);
end
eta.mean = nanmean(eta.activity,3);
eta.std = nanstd(eta.activity,[],3);

%%% now let's bootstrap if necessary
if yes_bootstrap == 1
    n_trials = size(eta.activity,3);
    eta.bootstrap = struct;
    eta.bootstrap.means = nan(size(eta.activity,1),size(eta.activity,2),bootstrap_n);        
    null_idx = [];
    for i = 1:size(eta_cell,1)
        sampleable = eta_cell{i,1}.bootstrap.sampleable;
        if size(sampleable,2)>size(sampleable,1)
            sampleable = sampleable';
        end
        null_idx = cat(1,null_idx,[ones(size(sampleable))*i sampleable]);
    end
    for i = 0:n_rounds-1
        % gather bootstrap_n/n_rounds * n_trials samples (with replacement)
        % we break it up to avoid memory problems
        bs_tmp_idx = datasample(null_idx,bootstrap_n*n_trials/n_rounds);         
        bs_data = nan(size(eta.activity,1),size(eta.activity,2),bootstrap_n*n_trials/n_rounds);
        unique_data = unique(bs_tmp_idx(:,1));
        % now slot that data in
        for d = 1:numel(unique_data)
            data_num = unique_data(d);
            data_idx = find(bs_tmp_idx(:,1)==unique_data(d));
            events_idx = bs_tmp_idx(data_idx,2);
            this_data = eta_cell{data_num,2};
            if ischar(this_data)
                this_data = load(this_data);
            end
            if isstruct(this_data)
                this_data = this_data.Fc; % just gonna assume it's an roi struct
            end
            if ~isempty(roi_nums)
                this_data = this_data(:,roi_nums);
            end
            if ~isempty(hp)
                this_data = highpass(this_data,hp(1),hp(2));
            end
            bs_eta = eventTriggeredAverage(this_data,events_idx,1,size(eta.activity,1),'nullDistr',0);
            bs_data(:,:,data_idx) = bs_eta.activity;
        end
        % now take the means of each chunk
        bs_idx_step = 1:n_trials:n_trials*bootstrap_n/n_rounds+1;
        for it = 1:bootstrap_n/n_rounds
            %disp([it it+(bootstrap_n*i/n_rounds)]);
            eta.bootstrap.means(:,:,it+(bootstrap_n*i/n_rounds)) = nanmean(bs_data(:,:,bs_idx_step(it):(bs_idx_step(it+1)-1)),3);
        end 
        %disp(i+1)
        fprintf('.')
    end
    fprintf('\n')

    %%% some stats
    eta.bootstrap.meansMu = nanmean(eta.bootstrap.means,3);
    eta.bootstrap.meansStd = nanstd(eta.bootstrap.means,[],3);
    eta.bootstrap.meansL = prctile(eta.bootstrap.means,2.5,3);
    eta.bootstrap.meansU = prctile(eta.bootstrap.means,97.5,3);
    eta.bootstrap.meansL1 = prctile(eta.bootstrap.means,0.5,3);
    eta.bootstrap.meansU1 = prctile(eta.bootstrap.means,99.5,3);
    eta.bootstrap.n = bootstrap_n;
    eta.bootstrap = rmfield(eta.bootstrap,'means');    
    % significant difference: -1 for sig neg, +1 for sig pos
    eta.bootstrap.meanSig = zeros(size(eta.mean));
    eta.bootstrap.meanSig(eta.mean<eta.bootstrap.meansL)=-1;
    eta.bootstrap.meanSig(eta.mean>eta.bootstrap.meansU)=1;
    % significant difference: -1 for sig neg, +1 for sig pos
    eta.bootstrap.meanSemSig = zeros(size(eta.mean));
    eta.bootstrap.meanSemSig(eta.mean+eta.std/sqrt(n_trials)<eta.bootstrap.meansL)=-1;
    eta.bootstrap.meanSemSig(eta.mean-eta.std/sqrt(n_trials)>eta.bootstrap.meansU)=1;        
end

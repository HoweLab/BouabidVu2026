% function output = get_reward_info(behav)
%
% This function takes as input a behavioral file 
% (the filepath, or the struct) and returns a struct with useful reward 
% information:
%   - rew_del:              the index of each reward delivery
%   - rew_size:             reward size (1 = sm, 2 = lg, 3 = md)
%   - rew_lick:             the index of the 1st lick after reward delivery
%   - time_since_rew_del:   the # s since the last reward delivery
%   - time_since_rew_con:   the # s since the last reward consumption lick
%
%
% updated 2/24/21 by Mai-Anh for general use 
% updated 6/1/21 by Mai-Anh for slightly more input flexibility
% updated 6/8/22 by Mai-Anh to add time_since_rew

function output = get_reward_info(behav)

% load if not a struct
if ~isstruct(behav)
    if ~endsWith(behav,'.mat')
        behav = [behav '.mat'];
    end
    behav = load(behav);
end

% delivery onsets
sr = round(1/mean(diff(behav.timestamp)));  % assume integer sampling rate
rew_on = vec(strfind(behav.reward',[0 1])+1); % delivery onsets
% if we have an unbinned datafile (not sampled to match neural recording), 
% assume rewards are at least 2s apart
if sr>1000 
    rew_keep = [1; diff(rew_on)>2000];
    rew_on = rew_on(rew_keep==1);
end

% reward size
if isfield(behav,'experimentSetup')
    rew_size = behav.experimentSetup.rew.rewSize;
elseif isfield(behav,'otherVars')
    rew_size = behav.otherVars.rew.rewSize;
elseif isfield(behav,'virmenInfo')
    if isfield(behav.virmenInfo.trials,'sizeSL')
        rew_size = behav.virmenInfo.trials.sizeSL;
    end
else 
    rew_size = [];
end

% consumption onsets: 1st lick after each reward 
lick_onset = nan(size(rew_on));
% all lick onsets
lick_on = find(behav.lick==1);
% consumption lick has to occur before the next reward to be considered
% part of current reward
rew_off = [rew_on(2:end)-1; numel(behav.reward)];
% reward-consumption activity
for r = 1:numel(rew_on)
    this_lick = lick_on(find(lick_on>(rew_on(r)),1,'first'));
    if ~isempty(this_lick) && this_lick<rew_off(r)
        lick_onset(r) = this_lick;
    end        
end

% time since last reward
time_since_rew_del = zeros(size(behav.timestamp));
for r = 1:numel(rew_on)
    this_rew = rew_on(r);
    time_since_rew_del(this_rew:end) = ...
        behav.timestamp(this_rew:end)-behav.timestamp(this_rew);
end

% time since last reward consumption
time_since_rew_con = zeros(size(behav.timestamp));
for r = 1:numel(lick_onset)
    this_rew = lick_onset(r);
    if ~isnan(this_rew)
        time_since_rew_con(this_rew:end) = ...
            behav.timestamp(this_rew:end)-behav.timestamp(this_rew);
    end
end

% compile output
output.rew_del = rew_on; % idx of rew delivery
output.rew_size = rew_size; % size
output.rew_lick = lick_onset; % idx of rew consumption (first lick after delivery)
output.time_since_rew_del = time_since_rew_del; % time(s) since most recent delivery
output.time_since_rew_con = time_since_rew_con; % time(s) since most recent consumption





% ETAstruct = eventTriggeredAverage(data,events,window1,window2,varargin);
%
% This function calculates the event-triggered activity and averages, as 
% well as a bootstrapped null distribution.
% The bootstrap null distribution is calculated by taking repeated samples
% of the same size as the events, and checking for significance against a
% percentile interval of the sample means.
%
% takes as input 
%   - data (n x p matrix): n = datapoints; p = diff timeseries (eg ROIs)
%   - events: a vector the indices of event occurrences
%   - window1: index of beginning of averaging window 
%       (<0 means before the event, >0 means after the event)
%   - window2: index of end of averaging window
%       (<0 means before the event, >0 means after the event)
%
% optional input
%   - otherExcl: other events to exclude from the bootstrap
%   - nullDistr: whether you want to estimate a bootstrapped null distrib 
%                   (0 = no, 1  = yes)
%   - bootstrapN: how many bootstrapped samples (default 5000)
%   - bootstrapSig: what %interval for significance (default .95)
%   - plotResults: whether you want to plot the results (default: 0 = no)
% 
%
% returns a struct (ETAstruct), with fields
%   - activity: event triggered activity: r x c x n
%       - r = another pt in the averaging window
%       - c = separate timeseries, eg, ROIs
%       - n = event occurrance 
%       - so if we gave it 2 ROIs, and want a window of 5 before and 
%         5 after, with 20 events happening, this matrix would be 11x2x20
%   - mean: the event triggered average: r x c
%   - std: event-triggered standard deviation: r x c
%   - events: the indices of the events included in the triggered average
%   - windowIdx: the relative indices of window
%   - bootstrap.means: the means from each of the bootstrapped samples 
%   - bootstrap.meansMu: the means of the sample means
%   - bootstrap.meansU: the upper of the 95th percentile (or whatever
%           percentile interval)
%   - bootstrap.meansL: the lower of the 95th percentile (or whatever
%           percentile interval)
%   - bootstrap.meanSig: a boolean indicating where the event-triggered average is
%              significant,i.e., mean falls outside the significance interval
%               you set of the values of the null distribution. For
%               example, if you set this interval to be 95%, then a value
%               of -1 means the triggered average mean falls below the
%               2.5th percentile and a value of 1 means that your triggered
%               average mean falls above the 97.5th percentile.
%   - bootstrap.meanSemSig: a boolean indicating where the event-triggered average is
%               significant,i.e., mean +/- SEM falls outside the significance interval
%               you set of the values of the null distribution. For
%               example, if you set this interval to be 95%, then a value
%               of -1 means the triggered average mean+sem falls below the
%               2.5th percentile and a value of 1 means that your triggered
%               average mean-sem falls above the 97.5th percentile.
%   - bootstrap.sampleable: the indices of parts of the timeseries that can
%               be sampled for the null distribution (i.e., not the
%               excluded timepoints)
%
% Mai-Anh Vu, 4/30/20
% updated 5/15/20
% updated 8/19/20 to check data size for single array data
% updated 1/21/21 by Mai-Anh to do better bootstrap
% updated 3/23/21 by Mai-Anh to allow option to plot output
% updated 3/25/21 by Mai-Anh: fixed commenting above
% updated 7/11/21 by Mai-Anh: option of saving out figures
% updated 9/27/22 by Mai-Anh: saves out the indices that are sampleable

function ETAstruct = eventTriggeredAverage(data,events,window1,window2,varargin)



%%%  parse optional inputs %%%
ip = inputParser;
ip.addParameter('otherExcl',[]);
ip.addParameter('nullDistr',1);
ip.addParameter('bootstrapN',5000);
ip.addParameter('bootstrapSig',.95);
ip.addParameter('plotResults',0);
ip.addParameter('saveDir',[]);
ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% some setup
if size(events,2)>size(events,1)
    events = events';
end
if size(data,1)==1 
    if size(data,2)>size(data,1)
        data = transpose(data);
    else
        disp('error: data is only 1 number.')
        return
    end
end
% output
ETAstruct = struct;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get event-related activity: index matrix strategy to avoid looping 
% the row indices (timepoints) of our data matrix
windowmat = repmat((window1:window2)',1,size(data,2),numel(events));
eventmat = permute(repmat(events,1,size(data,2),window2-window1+1),[3 2 1]);
% disp(size(windowmat))
% disp(size(eventmat))
rmat = windowmat+eventmat; 
rnan = rmat<1 | rmat>size(data,1); % indices falling outside our data
rmat(rnan) = 1;
% the column indices
cmat = repmat(1:size(data,2),window2-window1+1,1,numel(events));
% get linear indices
linmat = sub2ind(size(data),rmat,cmat);
% get the activity
ETAstruct.activity = data(linmat);
ETAstruct.activity(rnan) = nan;

% average
ETAstruct.mean = nanmean(ETAstruct.activity,3);
% standard deviation
ETAstruct.std = nanstd(ETAstruct.activity,[],3);

% the window indices
ETAstruct.events = events;
ETAstruct.windowIdx = transpose(window1:window2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bootstrapped null distribution

if nullDistr==1
    % some setup
    ETAstruct.bootstrap.means = nan(window2-window1+1,size(data,2),bootstrapN);
    % excluded events
    if ~isempty(otherExcl) && size(otherExcl,2)>size(otherExcl,1)
        otherExcl = otherExcl';
    end
    excludedEvents = unique([events; otherExcl]);
    includedEvents = (max([0 -window1])+1):(size(data,1)-window2);
    includedEvents = setdiff(includedEvents,excludedEvents);
    includedEvents = includedEvents(includedEvents>0);
    ETAstruct.bootstrap.sampleable = includedEvents;
    ETAstruct.bootstrap.sampleable = ETAstruct.bootstrap.sampleable(...
        ETAstruct.bootstrap.sampleable>(window2-window1+1));
    ETAstruct.bootstrap.sampleable = ETAstruct.bootstrap.sampleable(...
        ETAstruct.bootstrap.sampleable<(size(data,1)-window2+window1));
    

    for i = 1:bootstrapN
        % resample from our includedEvents with replacement
        bsevents = includedEvents(randi(numel(includedEvents),numel(events),1))';
        % the row indices (timepoints) of our data matrix
        windowmat = repmat((window1:window2)',1,size(data,2),numel(bsevents));
        eventmat = permute(repmat(bsevents,1,size(data,2),window2-window1+1),[3 2 1]);
        rmat = windowmat+eventmat; 
        % the column indices
        cmat = repmat(1:size(data,2),window2-window1+1,1,numel(bsevents));
        % get linear indices
        linmat = sub2ind(size(data),rmat,cmat);
        % activity
        temp = data(linmat);
        ETAstruct.bootstrap.means(:,:,i) = mean(temp,3);
    end
    % bootstrap average, , and upper and lower 95% percentiles
    ETAstruct.bootstrap.meansMu = nanmean(ETAstruct.bootstrap.means,3);
    ETAstruct.bootstrap.meansU = quantile(ETAstruct.bootstrap.means,1-(1-bootstrapSig)/2,3);
    ETAstruct.bootstrap.meansL = quantile(ETAstruct.bootstrap.means,(1-bootstrapSig)/2,3);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % significant difference: -1 for sig neg, +1 for sig pos
    ETAstruct.bootstrap.meanSig = zeros(size(ETAstruct.mean));
    ETAstruct.bootstrap.meanSig(ETAstruct.mean<ETAstruct.bootstrap.meansL)=-1;
    ETAstruct.bootstrap.meanSig(ETAstruct.mean>ETAstruct.bootstrap.meansU)=1;
    % significant difference: -1 for sig neg, +1 for sig pos
    ETAstruct.bootstrap.meanSemSig = zeros(size(ETAstruct.mean));
    ETAstruct.bootstrap.meanSemSig(ETAstruct.mean+ETAstruct.std/sqrt(numel(events))<ETAstruct.bootstrap.meansL)=-1;
    ETAstruct.bootstrap.meanSemSig(ETAstruct.mean-ETAstruct.std/sqrt(numel(events))>ETAstruct.bootstrap.meansU)=1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot output
if plotResults == 1
    % get subplot arrangement
    nData = size(data,2);
    [nRows, nCols, ~] = nSubplots(nData);
    
    % loop over each data column, starting a new figure when necessary
    k = 0;
    fignum = 0;
    for d = 1:nData
        % open a new (full-sized) figure if necessary
        if k==0 
            figure('units','normalized','outerposition',[0 0 1 1],'visible','on')
            fignum = fignum+1; % increment fig num
        end
        k=k+1; % increment k
        
        % now plot
        subplot(nRows,nCols,k)
        hold on
        x = ETAstruct.windowIdx;
        
        % if there's bootstrap info
        if nullDistr==1
            
            % color
            thisColor = [.7 .7 .7];
            
            % mean, and upper (U) and lower (L) bootstrap null distr int
            mu = ETAstruct.bootstrap.meansMu(:,d);
            muU = ETAstruct.bootstrap.meansU(:,d);
            muL = ETAstruct.bootstrap.meansL(:,d);

            % plot sem then mean
            fill([x; flipud(x)],[muU; flipud(muL)],thisColor,'FaceAlpha',.3,'LineStyle','none')
            plot(x,mu,'-','Color',thisColor);
        end
        
        % plot the event triggered SEM & mean
        thisColor = lines(1);
        sigColor = lines(7);
        sigColor = sigColor(7,:);
        mu = ETAstruct.mean(:,d);
        sem = ETAstruct.std(:,d)/sqrt(size(ETAstruct.activity,3));
        % plot sem and mean
        fill([x; flipud(x)],[mu-sem; flipud(mu+sem)],thisColor,'FaceAlpha',.3,'LineStyle','none')
        plot(x,mu,'-','Color',thisColor);
        
        % if there's bootstrap info, indicate significance
        if nullDistr==1
            % flag significant points with red dot
            sig = ETAstruct.bootstrap.meanSemSig(:,d);
            plot(x(sig~=0),mu(sig~=0),'.','Color',sigColor)
        end
        
        % informative lines and text
        % put on a x=0 line if it's useful
        ys = get(gca,'YLim');
        if ismember(0,x)
            plot([0 0],ys,'-k');
        end
        % put on a y=0 line if it's useful
        if ys(1)<0 %& ys(2)>0
            plot([x(1) x(end)],[0 0],'-k')
        end        
        % axes & titles
        set(gca,'XLim',[x(1) x(end)],'YLim',ys)
        xlabel('idx')
        title(['Data ' num2str(d)]);       
        
        % if we're at the end of the figure (max #subplots, reset k to next
        % loop will open a new fig)
        if k==nRows*nCols || d == nData
            k = 0;
            if ~isempty(saveDir)
                saveas(gcf,fullfile(saveDir,['fig' sprintf('%02d',fignum)]),'fig');
            end
        end 
    end   
end
end

%%%%%%%%%%%%%%%%%%%%%%%
% subplot configuration
function [subplotRows, subplotCols, nFigs] = nSubplots(n)

    % some posssible subplot arrangments (>20 will open multiple figs)
    arr = [
        1 1;...
        2 2;...
        2 3;...
        3 3;...
        3 4;...
        4 4;...
        4 5];
    arr_n = prod(arr,2);
    which_arr = find(arr_n >= n,1,'first');
    if isempty(which_arr)
        which_arr = size(arr,1);
    end
    nFigs = ceil(n/arr_n(which_arr));

    subplotRows = arr(which_arr,1);
    subplotCols = arr(which_arr,2);

end


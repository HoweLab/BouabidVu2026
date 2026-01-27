% [output,retry] = split_and_bin_grid_behav_files(behavFilepath,varargin)
%
% This function accomplishes 2 things, based on TTLs sent on
% each frame acquistion from an imaging device, such as Hamamatsu or
% Scanbox.
% 1. It splits them into movies, based on gaps between bouts of TTL inputs.
% 2. It and bins the data according to these TTLs (e.g., aligns
% them to all the received TTLs from the camera). 
%
% More specifically: this function reads in the path to a struct of data
% (such as one acquired via NIDAQ and Mai-Anh's GUI.m), with TTLs coming,
% once per frame,from an imaging source denoted by fields 'ttlIn1' and/or 
% 'ttlIn2'. If there isn't any ttl input but you want to align the two as if
% there are, then set the optional input 'noNeural' to 1 (and read details
% below). It then saves out all non-numeric data as is. Numeric data are 
% handled as outlined below. Virmen data is also treated the same way.
%
% Definition of bins by TTLs:
% By default, this centers each "bin" in the middle of each frame (i.e., a
% bin starts at the rising edge of a TTL and ends before the next TTL 
% rising edge. We use 'VSYNC' for our output trigger option, which means
% that the TTL gets sent 500us (what we set as our delay typically) after 
% the first line is finished exposing and has begin to be written out. So 
% this means are bins begin as soon as the first line is exposed. 
% eg: split_and_bin_grid_behav_files(behavFilepath,'centering','middle');
% The other centering option is "beginning", where the bin is centered on
% the beginning of the frame (i.e., the bin is centered around the TTL
% rising edge, such that the bin begins halfway between the previous TTL, 
% and ends halfway between the current TTL and next TTL).
% eg: split_and_bin_grid_behav_files(behavFilepath,'centering','beginning');
%
% Notes on how variables are binned: 
% Binary variables are handled in two ways: binary and count. 
% When handled as binary variables, we assign a 1 if any value within the 
% bin was a 1 (i.e., we take the maximum). When handled as count variables, 
% denoted in the output with '_count', we count the number of 1s in the bin,
% (i.e., we take the sum). 
% Continuous variables are meaned.
%
% (5/12/25) The spacing between the TTLs will be checked, and if further 
% attention is warranted, the check_ttl_spacing spacing field will be set 
% to TRUE, and otherwise FALSE. It will also display a message to the 
% command line.
% 
% Required input:
% the path to the behavioral file as a string
%
% Optional input:
% 'centering' -- set to 'beginning' or 'middle' to center the bins in the
%       beginning of the frame or the center of it (see description above). 
%       Default is "middle".
% 'splitThresh' -- when splitting one long acquisition into separate
%       movies, the amt of time to use as a detector of "gaps" in between
%       acquisitions/focusing. Default is 250 (ms).
% 'virmenStartTTL' -- which ttlOut (Nidaq) to use to align with the onset 
%       of the first trial of the ViRMEn task. Note that a TTL is sent at 
%       the start of each new trial (phase=1). Default is 3, for the cue 
%        rotation task. The first 2 TTLs sent are to signal the start of 
%       the acquisition. To figure this out, look for the double TTL. The
%       virmenStartTTL will be the next one.
% 'noNeural' -- if there is no neural data, set this to 1 (default is 0).
%       If this is set to 1, then we'll just "make up" TTLs to align
%       everything to
% 'noNeural_Hz' -- in the case that noNeural (above) is set to 1, what
%       sampling frequency would you like to align to? Default is 50Hz.
% 'ttlInField' -- in the case of splitting it by a specific ttlIn field,
%       provide that here, e.g., {'ttlIn1'}, or {'ttlIn1','ttlIn2'}.
%       default is {'ttlIn1','ttlIn2'}
%
% Output:
% This will output .mat file(s) to the same location as the original
% behavioral file and will append "movie#" to the end, depending on
% how many "movies", i.e., bouts of TTLs, are detected. This function will
% also output a figure showing where movie boundaries were detected.
%
% 
%
% Mai-Anh Vu
% updated 11/27/2019
% updated 5/14/2020
% updated 8/14/2020 to handle TTLs of 2 different orientations
% updated 8/18/2020 to save out binned binary variables both as binary, and as count variables
% updated 8/22/2020 documentation above
% updated 12/14/2020 
% updated 1/26/2021 to allow alignment of behavioral and virmen data without
% neural data by just putting on TTLs
% updated 6/3/2021 to handle TTLs that go all the way to end of recording
% updated 5/12/2025 to change how we check TTL spacing
%
% obsolete previous updates:
% updated 5/21/21 to detect and ignore erroneous TTLs
% updated 9/21/21 to allow optional input determining the erroneous TTL
%   threshold (percent difference from mode spacing)
% updated 8/22/22 to allow turning off of erroneous TTL detection (this
%   will automatically turn off if noNeural is set to 1)
% updated 10/06/22 (Zack) last commits convert ttl_ignored from binary array
% to index array (change it back in this commits)


function [output, retry] = split_and_bin_grid_behav_files(behavFilepath,varargin)
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  parse inputs %%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    ip = inputParser;
    ip.addParameter('centering','middle');
    ip.addParameter('splitThresh',250);
    ip.addParameter('virmenStartTTL',3);
    ip.addParameter('noNeural',0);
    ip.addParameter('noNeural_Hz',50);
    ip.addParameter('ttlInField',{'ttlIn1','ttlIn2'});
    ip.parse(varargin{:});
    for j=fields(ip.Results)'
        eval([j{1} '=ip.Results.' j{1} ';']);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  load data % setup %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~endsWith(behavFilepath,'.mat')
        behavFilepath = [behavFilepath '.mat'];
    end
    behav = load(behavFilepath);
    if strncmp(centering,'beginning',3)
        centered = 0;
    elseif strncmp(centering,'middle',3)
        centered = 1;
    else
        disp('  NOTE: The variable ''centering'' was not set to either')
        disp('   ''middle'' or ''beginning''. This script will default ')
        disp('   to ''middle''. Type "help split_and_bin_grid_behav_files"');
        disp('   for more info.');

        centered = 1;
    end
    
    % ttlIn fields: get rid of the field if it has less than 1% of the ttls
    % from the other field
    sumField = zeros(size(ttlInField));
    for ti = 1:numel(ttlInField)
        tiSum(ti) = sum(behav.(ttlInField{ti}));
    end
    ttlInField(tiSum<max(tiSum)/100) = [];
    retry = zeros(size(ttlInField));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  parse apart movies %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % if we did record neural data, parse apart the movies
    if noNeural==0
        % note: there are multiple movies per acquisition, so let's parse them
        % apart
        % frame TTLs (look for the upticks in the ttlIn trace)
        if isfield(behav,'ttlIn_Scanbox') && ~isfield(behav,'ttlIn2')
            behav.ttlIn2 = behav.ttlIn_Scanbox;
        end
        if isfield(behav,'ttlIn_Hamamatsu') && ~isfield(behav,'ttlIn1')
            behav.ttlIn1 = behav.ttlIn_Hamamatsu;
        end
        potentialTTLinFields = ttlInField;
        keepFields = zeros(size(potentialTTLinFields));
        for f = 1:numel(potentialTTLinFields)
            if isfield(behav,potentialTTLinFields{f})
                if sum(behav.(potentialTTLinFields{f}))>0
                    keepFields(f) = 1;
                end        
            end
        end    
        if sum(keepFields)==0
            disp('  NOTE: Sorry, there seems to be no TTL input.')
            return
        end
        whichFields = find(keepFields==1);
    else
        potentialTTLinFields = {'ttlIn000'};
        whichFields = 1;
        behav.ttlIn000 = zeros(size(behav.ttlOut));
        % make up some TTLs, spaced according to noNeural_Hz
        % leaving 1s of buffer on either end
        behav.ttlIn000(2000:round((2000/noNeural_Hz)):end-2000) = 1; 
    end
    % loop over input TTL fields
    for ti = 1:numel(whichFields)
        ttlInFieldName = potentialTTLinFields{whichFields(ti)};
        ttlField = behav.(ttlInFieldName); 
        
        %%%%% SPLIT INTO MOVIES
        % let's first figure out if TTL orientation is pos or neg
        % if % more 0s than 1s -> pos -- look for rising edge
        if sum(ttlField==0)>sum(ttlField==1) 
            ttlInOn = vec(strfind(transpose(vec(ttlField)),[0 1])+1);
        else % more 1s than 0s -> neg -- look for falling edge
            ttlInOn = vec(strfind(transpose(vec(ttlField)),[1 0])+1);
        end
        % let's separate out movies
        ttlSpaces = [1; diff(ttlInOn)];
        % look for gaps to indicate starts/ends of recording
        fs = mean(diff(behav.timestamp)); % behavioral file sampling frequency
        ttlBigGaps = find(ttlSpaces>(splitThresh/1000/fs));
        if ~isempty(ttlBigGaps)
            movieStarts = [ttlInOn(1); ttlInOn(ttlBigGaps)]; 
            movieEnds = [ttlInOn(ttlBigGaps-1); ttlInOn(end)]; 
        else
            movieStarts = ttlInOn(1);
            movieEnds = ttlInOn(end);
        end
        ignore = movieEnds-movieStarts==0;
        movieEnds = movieEnds(~ignore);
        movieStarts = movieStarts(~ignore);
    
        % check movie identification
        if noNeural==0
            figure
            plot(ttlField)
            hold on
            for m = 1:numel(movieStarts)
                plot([movieEnds(m) movieEnds(m)],[0 1],'r','LineWidth',2)
                plot([movieStarts(m) movieStarts(m)],[0 1],'g','LineWidth',2)   
            end
            title('movie boundaries')
            savefig(gcf,[behavFilepath(1:end-4) '_' ttlInFieldName '_movieBoundaries.fig'])
            hold off
            close all
        end
    
        % number of frames per movie
        movieLengths = nan(numel(movieStarts),1);
        starti = nan(numel(movieStarts),1);
        endi = nan(numel(movieStarts),1);
        % put a buffer of 50ms (100 datapts)before start and after end 
        % (or less, if we reach the beginning or end of the recording)
        buf = 100; 
        for m = 1:numel(movieStarts)  
            starti(m) = max([movieStarts(m)-buf 1]);
            endi(m) = min([movieEnds(m)+buf numel(ttlField)]);
            temp = ttlField(starti(m):endi(m));
            movieLengths(m) = sum(diff(temp)==1);
        end
        
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% output setup %%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % separate data fields from information fields
        structfields = fieldnames(behav);
        infoFields = structfields(structfun(@(x) ~isnumeric(x),behav));
        dataFields = structfields(structfun(@(x) isnumeric(x),behav));
        % if there's virmen data
        if isfield(behav,'virmen')
            infoFields = infoFields(~strncmp(infoFields,'virmen',6));
            virmenDataCount = max(structfun(@(x) numel(x),behav.virmen));
            virmenFields = fieldnames(behav.virmen);
            virmenInfoFields = virmenFields(structfun(@(x) numel(x)<virmenDataCount, behav.virmen));
            virmenDataFields = virmenFields(structfun(@(x) numel(x)==virmenDataCount, behav.virmen));
        end
    
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% fill in output %%%%%
        %%%%%%%% and save %%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % loop over movies
        for m = 1:numel(movieStarts)    
            % preallocate output
            output = struct;
            % fill in the informational fields
            for f = 1:numel(infoFields)
                output.(infoFields{f}) = behav.(infoFields{f});
            end
            % preallocate numeric fields based on #frames for each movie  
            for f = 1:numel(dataFields)
                output.(dataFields{f}) = nan(movieLengths(m),1);
            end
            % now let's find where the frames are  
            ttlIn = ttlField(starti(m):endi(m));
            % find each uptick for each frame
            frameOn = find(diff(ttlIn)==1)+1;
            % av bin size
            avBin = round(mean(diff(frameOn)));    
            if centered == 1 % TTLs start the bins
                binBounds = frameOn;
                binBounds = [binBounds; binBounds(end)+avBin];        
            else % TTLs are in the middle of the bins 
                % find the midpoints between frames
                mdpts = round(mean([frameOn(1:end-1) frameOn(2:end)],2));
                binBounds = [frameOn(1)-round(avBin/2); mdpts; frameOn(end)+round(avBin/2)];
            end
            binBounds = binBounds(binBounds<=numel(ttlIn));
            binBounds = unique(binBounds); % in case there are any duplicates
    
            % loop over fields
            dataFields_without_timestamp = dataFields(~contains(dataFields,"timestamp"));
            % set up timestamp by taking exact ttl onset
            this_timestamp = behav.timestamp(starti(m):endi(m));
            for t = 1:numel(binBounds)-1              
                output.timestamp(t) = this_timestamp(binBounds(t));
            end
            % set up rest non timestamp fields
            for f = 1:numel(dataFields_without_timestamp)
                % load this field (with start/end buffer)
                thisField = behav.(dataFields_without_timestamp{f})(starti(m):endi(m));
                % check if it's a binary variable
                isbin = sum(double(logical(thisField))-thisField)==0;
                % loop over frames
                for t = 1:numel(binBounds)-1              
                    if isbin % for binary variables, if there was a one in that bin, assign 1
                        output.(dataFields_without_timestamp{f})(t) = nanmax(thisField(binBounds(t):(binBounds(t+1)-1)));
                        output.([dataFields_without_timestamp{f} '_count'])(t,1) = nansum(thisField(binBounds(t):binBounds(t+1)-1));
                    else % for continuous variables, take the mean
                        output.(dataFields_without_timestamp{f})(t) = nanmean(thisField(binBounds(t):binBounds(t+1)-1));
                    end
                end
            end
            
            %%%% check for suspicious TTL spacing
            ttlInFieldName = potentialTTLinFields{whichFields(ti)};
            ttl_check = check_TTL_spacing(output,'ttl_channels',{ttlInFieldName});
            if ttl_check.(ttlInFieldName) == 1 % we got flagged to check this
                output.check_ttl_spacing = true;
                disp('  ALERT: you may need to double check the timing of')
                disp(['   the TTLs in ' ttlInFieldName  '. See:'])
                disp('    - check_TTL_spacing.m')
                disp('    - handle_erroneous_TTLs.m')
                disp('   This has also been noted in the output field') 
                disp('   check_ttl_spacing, which has been set to TRUE.')
            else
                output.check_ttl_spacing = false;
            end
            %%%%%%%%%
            % now let's handle virmen data, if applicable
            if isfield(behav,'virmen')
                behavTTLoffset = find(diff(behav.ttlOut)==1)+1;
                % use the second one: the first is the initial TTL
                behavTTLoffset = behav.timestamp(behavTTLoffset(virmenStartTTL)); 
                virmenTTLoffset = behav.virmen.phase==1;
                virmenTTLoffset = find(diff(virmenTTLoffset)==1)+1;
                virmenTTLoffset = behav.virmen.timestamp(virmenTTLoffset(1));
    
                output.virmenInfo = struct;
                % fill in information fields
                for f = 1:numel(virmenInfoFields)
                    output.virmenInfo.(virmenInfoFields{f}) = behav.virmen.(virmenInfoFields{f});
                end
                 % preallocate numeric fields based on #frames for each movie  
                for f = 1:numel(virmenDataFields)
                    output.(['virmen_' (virmenDataFields{f})]) = nan(movieLengths(m),1);
                end
                % let's get the timestamp (s) of each frame
                %zero-ing out at the first TTL out (which lines up with virmen
    
                timestamps = behav.timestamp(starti(m):endi(m))-behavTTLoffset;
                ttlTimestamps = timestamps(frameOn);        
                % av bin size
                avTimebin = mean(diff(ttlTimestamps));            
                if centered == 1 % if we're centering the bins in the middle of the frame
                    timebinBounds = ttlTimestamps;
                    timebinBounds = [timebinBounds; timebinBounds(end)+avTimebin];
                else % or if we're centering the bins on the beginning of the frame
                    % find the midpoints between frames
                    timeMdpts = mean([ttlTimestamps(1:end-1) ttlTimestamps(2:end)],2);            
                    timebinBounds = [ttlTimestamps(1)-avTimebin/2; timeMdpts; ttlTimestamps(end)+round(avTimebin/2)];
                end    
                % now get the indices of these bounds
                virmen_starti = find((behav.virmen.timestamp-virmenTTLoffset)>=(behav.timestamp(starti(m))-behavTTLoffset),1,'first');
                virmen_endi = find((behav.virmen.timestamp-virmenTTLoffset)<=(behav.timestamp(endi(m))-behavTTLoffset),1,'last');
                virmen_binBounds = nan(numel(timebinBounds),1);
                virmen_timestamp = behav.virmen.timestamp(virmen_starti:virmen_endi)-virmenTTLoffset;       
                for b = 1:numel(timebinBounds)
                    tmp = find(virmen_timestamp>=timebinBounds(b),1,'first');
                    if ~isempty(tmp)
                        virmen_binBounds(b) = tmp;
                    else
                        virmen_binBounds(b) = numel(virmen_timestamp);
                    end
                end            
                % loop over virmen data fields
                for f = 1:numel(virmenDataFields)
                    % load this field, and crop it down to the chunk of time of
                    % interest, with buffer            
                    wholeField = behav.virmen.(virmenDataFields{f});
                    thisField = behav.virmen.(virmenDataFields{f})(virmen_starti:virmen_endi);
                    % check if it's a binary variable
                    isbin = sum(double(logical(wholeField))-wholeField)==0;
                    % loop over frames
                    for t = 1:numel(frameOn)
                        tmp = thisField(virmen_binBounds(t):virmen_binBounds(t+1)-1);
                        if ~isempty(tmp)
                            if isbin % for binary variables, if there was a one in that bin, assign 1
                                output.(['virmen_' (virmenDataFields{f})])(t) = nanmax(tmp);
                                output.(['virmen_' (virmenDataFields{f}) '_count'])(t) = nansum(tmp);
                            else % for continuous variables, take the mean
                                output.(['virmen_' (virmenDataFields{f})])(t) = nanmean(tmp);
                            end
                        end
                    end
                end
                % let's do a check, comparing reward signal and virmen phase 3
                % reward sig
                rew = [0; diff(output.reward>0)];
                rew = rew==1;
                rt = output.timestamp(rew==1);
                % virmen phase 3
                phase3 = ceil(output.virmen_phase)==3;
                phase3 = [0;diff(phase3)];
                phase3 = phase3==1;
                vt = output.timestamp(phase3==1);
                % now check mean diff in timestamp
                if numel(vt)<numel(rt)
                    dt = abs(vt-repmat(rt',numel(vt),1));
                else
                    dt = abs(rt-repmat(vt',numel(rt),1));
                end
                dt_mean = mean(min(dt,[],2));
                if dt_mean>0.1 % if, on average, we're off by more than 100ms, flag
                    retry(ti) = 1;
                    disp('  *** ALERT: output saved, but misalignment likely.')
                    disp('   Check the alignment, and if needed, try changing')
                    disp('   ''virmenStartTTL'' and run again.')
                    disp('   Suggestion: first try increasing it by 1.')
                end
            end
            outFilepath = [behavFilepath(1:end-4) '_' ttlInFieldName '_movie' num2str(m) '.mat'];
            save(outFilepath,'-struct','output')       
            
        end
    end
end


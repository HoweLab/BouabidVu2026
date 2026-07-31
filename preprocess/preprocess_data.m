%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preprocess_data
%%%%%%%%%%%%%%%%%%%%%%%%%%%

function output = preprocess_data(path_roi1,path_roi2,path_behav1,path_behav2,varargin)
    % takes as input paths to corresponding roi files and behav files from
    % a single recording. if only a 1-channel recording, leave path_roi2
    % and path_behav2 empty ([]).
    
    %%%  parse optional inputs %%%
    ip = inputParser;
    ip.addParameter('light_blink_artifact',[0 0]); % artifacts from IR light blinking due to overheating    
    ip.addParameter('channel_names',{'roi1','roi2'}); % same order as input
    ip.addParameter('hp_hz',[0.3 0.1]); % high pass filter freq: default ch1 = 470nm = 0.3Hz, ch2 = 570nm = 0.1Hz
    ip.parse(varargin{:});
    for j=fields(ip.Results)'
        eval([j{1} '=ip.Results.' j{1} ';']);
    end
    
    
    % 1. align 2 channels if necessary
    if isempty(path_roi2) || isempty(path_behav2)
        output = struct;
        output.(channel_names{1}) = load(path_roi1);
        output.(['behav_' channel_names{1}]) = load(path_behav1);
    else
        output = align_2ch(path_roi1,path_roi2,path_behav1,path_behav2,...
            'channel_names',channel_names);        
    end        
    
    % roi deltaF/F preprocessing    
    for w = 1:sum(startsWith(fieldnames(output),'behav'))
        roiF = output.(channel_names{w}).F;
        sr = round(1/nanmean(diff(output.(['behav_' channel_names{w}]).timestamp)));
        
%         % sliding 8th percentile baseline
%         [Fc,F_baseline,~] = FtoFc(roiF,540);            
%         output.(channel_names{w}).FtoFcWindow = sr*30;
%         output.(channel_names{w}).Fc = Fc;
%         output.(channel_names{w}).F_baseline = F_baseline;

        % 1. calculate DFF from 2-term exponential baseline (capture fast
        % initial bleaching decay, and then slower decay)
        [Fc,F_baseline] = FtoFc_exp(roiF,'exp_model','exp2');
        output.(channel_names{w}).F_baseline_exp = F_baseline;  
        output.(channel_names{w}).Fc_exp = Fc;                

        % 2. highpass filter
        hp_cutoff = hp_hz(w);
        hp_steepness = 0.8;
        try
            output.(channel_names{w}).Fc_exp_hp = highpass(output.(channel_names{w}).Fc_exp,...
                hp_cutoff,sr,'ImpulseResponse','fir','steepness',hp_steepness);
        catch exception
            disp('     could not high-pass filter')
            output.(channel_names{w}).Fc_exp_hp = [];
        end  
        output.(channel_names{w}).hp_steepness = hp_steepness;
        output.(channel_names{w}).hp_hz = hp_hz(w);
                
    end
    
    % 3. artifact detection (if both channels are present; uses information from both channels)
    if isfield(output,channel_names{1}) && isfield(output,channel_names{2})
        art_results = classifyMultiChannelArtifacts(...
            output.(channel_names{1}).Fc_exp,output.(channel_names{2}).Fc_exp,...
            'Type1ZThresh',15,'Type1MinSlope',8,...
            'ChannelNames',channel_names,'Type2Likely',light_blink_artifact);
        for w = 1:numel(channel_names)
            output.(channel_names{w}).artifact_mask = art_results.(channel_names{w}).autoRemovedMask;
        end
    end
    
    % 4. whether or not that day has signal, based on # transients: 
    % use highpass filtered data (with artifact masked out) for this
    for w = 1:sum(startsWith(fieldnames(output),'behav'))
        if ~isempty(output.(channel_names{w}).Fc_exp_hp)
            this_sig_hp = output.(channel_names{w}).Fc_exp_hp;
            this_sig_hp(output.(channel_names{w}).artifact_mask) = nan;
            output.(channel_names{w}).sig = get_signal_tr(this_sig_hp);
        else
            output.(channel_names{w}).sig = 0;
        end        
    end
    
end
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% align_2ch
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output = align_2ch(path_roi1,path_roi2,path_behav1,path_behav2,varargin)
    % Take as input paths of roi data files corresponding to ttlIn1 and 
    % ttlIn2 NIDAQ inputs to behavior file (e.g., 470nm, 570nm), behav 
    % files already aligned to TTLs (ttl1, ttl2). Aligns them based on 
    % timestamp, truncate and ignore non-simultaneous frames.
    %
    % Output a struct
    %
    % modified from code by Liangzhu Zhang
    
    %%%  parse optional inputs %%%
    ip = inputParser;
    ip.addParameter('channel_names',{'roi1','roi2'}); % same order as input
    ip.addParameter('hp_hz',[0.3 0.1]); % high pass filter freq: default ch1 = 470nm = 0.3Hz, ch2 = 570nm = 0.1Hz
    ip.parse(varargin{:});
    for j=fields(ip.Results)'
        eval([j{1} '=ip.Results.' j{1} ';']);
    end

    roi1 = load(path_roi1); 
    roi2 = load(path_roi2); 
    behav1 = load(path_behav1); 
    behav2 = load(path_behav2);
    % just find common timestamps and align em...
    stamp1 = behav1.timestamp(~isnan(behav1.timestamp));
    stamp2 = behav2.timestamp(~isnan(behav2.timestamp));
    
    stamp1_se = stamp1([1,end]);
    stamp2_se = stamp2([1,end]);
    
    truncate_se = [max(stamp1_se(1),stamp2_se(1)),min(stamp1_se(2),stamp2_se(2))];
    
    [~,I1start] = min(abs(behav1.timestamp-truncate_se(1)));
    [~,I1end] = min(abs(behav1.timestamp-truncate_se(2)));
    [~,I2start] = min(abs(behav2.timestamp-truncate_se(1)));
    [~,I2end] = min(abs(behav2.timestamp-truncate_se(2)));

    len_F1 = size(roi1.F,1);
    len_F2 = size(roi2.F,1);
    I_length = min([I1end-I1start,I2end-I2start,len_F1-I1start,len_F2-I2start]);
    I1end = I1start + I_length;
    I2end = I2start + I_length;    
    
    % slicing F-related fields
    F_fields = fieldnames(roi1);
    F_fields = F_fields(structfun(@(x) size(x,1) == size(roi1.F,1),roi1));
    for f = 1:numel(F_fields)
        roi1.(F_fields{f}) = roi1.(F_fields{f})(I1start:I1end,:);
    end    
    F_fields = fieldnames(roi2);
    F_fields = F_fields(structfun(@(x) size(x,1) == size(roi2.F,1),roi2));
    for f = 1:numel(F_fields)
        roi2.(F_fields{f}) = roi2.(F_fields{f})(I2start:I2end,:);
    end
    
    % slicing behav fields
    try        
        n_frames_2 = size(behav2.timestamp,1);
        for field = string(fields(behav2)')
            if isnumeric(behav2.(field)) && size(behav2.(field),1)==n_frames_2
                behav2.(field) = behav2.(field)(I2start:I2end);
            elseif isnumeric(behav2.(field)) && size(behav2.(field),1)==n_frames_2-1
                behav2.(field) = behav2.(field)(I2start:I2end-1);
            end
        end
        n_frames_1 = size(behav1.timestamp,1);
        for field = string(fields(behav1)')
            if isnumeric(behav1.(field)) && size(behav1.(field),1)==n_frames_1
                behav1.(field) = behav1.(field)(I1start:I1end);
            elseif isnumeric(behav1.(field)) && size(behav1.(field),1)==n_frames_1-1
                behav1.(field) = behav1.(field)(I1start:I1end-1);
            end
        end
    catch err
        error("Error Probably due to weird fields in behav files, check it.\n%s -> %s\n",err.identifier,err.message)
    end

    % output
    output.(channel_names{1}) = roi1;
    output.(channel_names{2}) = roi2;
    output.(['behav_' channel_names{1}]) = behav1;
    output.(['behav_' channel_names{2}]) = behav2;    
    output.([channel_names{1} '_idx']) = [I1start I1end];
    output.([channel_names{2} '_idx']) = [I2start I2end];
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FtoFc
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ Fc, scale, center ] = FtoFc( F, scale_window )
    %FTOFC Normalizes by the 8th percentile in a sliding window and subtracts
    %   the median to create the modified DFoF

    if ~exist('scale_window','var'), scale_window = 1e3; end

    scale = zeros(size(F));
    for i=1:size(scale,1)
        scale(i,:) = quantile(F(max(i-scale_window,1):min(i+scale_window,size(F,1)),:),0.08);
    end
    Fc = F./scale;
    center = median(Fc);
    Fc = bsxfun(@minus,Fc,center);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FtoFc_exp
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ Fc, scale, center ] = FtoFc_exp(F,varargin)
    % An alternative to the FtoFc function to calculate DFF.
    % FtoFc.m calculates the baseline (B) using a sliding 8th percentile 
    % window, and then normalizes F to that B in the following way:
    % (F./B) - median(F./B)
    %
    % This function calculates the baseline B by fitting 2-term exponential
    % function to the F instead.
    %
    % Mai-Anh Vu
    % 5/16/2023
    % edited 12/17/2024 so you can feed a different exponential model using
    % input 'exp_model'

    %%%  parse optional inputs %%%
    ip = inputParser;
    ip.addParameter('exp_model','exp2'); % other possibility is exp1
    ip.parse(varargin{:});
    for j=fields(ip.Results)'
        eval([j{1} '=ip.Results.' j{1} ';']);
    end

    % now fit the exponential to calculate the baseline
    y_all = F;
    scale = nan(size(y_all));
    x = 1:size(F,1);
    x = x(:);
    for i = 1:size(y_all,2)
        this_y = y_all(:,i);
        this_x = x(:);
        are_nan = isnan(this_y) | isnan(this_x);  
        this_x = this_x(~are_nan);
        this_y = this_y(~are_nan);    
        mdl = fit(this_x,this_y,exp_model);       
        this_scale = mdl(this_x);
        scale(~are_nan,i) = this_scale;
    end

    % DFF normalization
    Fc = F./scale;
    center = nanmedian(Fc);
    Fc = Fc - center;
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


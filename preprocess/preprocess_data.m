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
    ip.addParameter('z_motion_artifact',0); % usually animal-specific; whether z-motion is likely
    ip.addParameter('channel_names',{'roi1','roi2'}); % same order as input
    ip.addParameter('hp_hz',[0.3 0.1]); % high pass filter freq: default ch1 = 470nm = 0.3Hz, ch2 = 570nm = 0.1Hz    
    ip.addParameter('Type1ZThresh',10); % amplitude z threshold for identifying error type 1
    ip.addParameter('Type1MinSlope',10); % slope z threshold for identifying error type 1
    ip.addParameter('Type3MinDurSRFrac',2); % sr/Type3MinDurationFrac is the min duration threshold for error type 3
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
    
    % roi deltaF/F preprocessing: initial
    for c = 1:sum(startsWith(fieldnames(output),'behav'))
        roiF = output.(channel_names{c}).F;
        sr = round(1/nanmean(diff(output.(['behav_' channel_names{c}]).timestamp)));
        
        % 1. DFF: sliding 8th percentile baseline
        sr_window = sr*30/2; % 30s window (+/- 15s)
        [Fc,F_baseline,~] = FtoFc(roiF,'F_to_baseline',despike(roiF),'scale_window',sr_window);            
        output.(channel_names{c}).FtoFcWindow = sr_window;  
        output.(channel_names{c}).Fc = Fc;
        output.(channel_names{c}).F_baseline = F_baseline;       
                
    end
    
    % 2. artifact detection (if both channels are present; uses information from both channels)
    if isfield(output,channel_names{1}) && isfield(output,channel_names{2})
        % 30s sliding 8th percentile window
        art_results = classifyMultiChannelArtifacts(...
            output.(channel_names{1}).Fc,output.(channel_names{2}).Fc,...
            'ChannelNames',channel_names,...
            'Type1ZThresh',Type1ZThresh,'Type1MinSlope',Type1MinSlope,...   % Type1 artifact settings
            'Type2Likely',light_blink_artifact,...                          % Type2 artifact settings
            'Type3MinDuration',sr/Type3MinDurSRFrac,...                     % Type3 artifact settings
            'Type3Likely',z_motion_artifact);    
        for c = 1:numel(channel_names)
            for a = 1:3
                output.(channel_names{c}).(['artifact' num2str(a)]) = ...
                    art_results.(channel_names{c}).(['type' num2str(a) 'Mask']);
            end
            output.(channel_names{c}).artifact_mask = art_results.(channel_names{c}).autoRemovedMask;
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
% despike
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F_despiked = despike(F,varargin)
    
    %%%  parse optional inputs %%%
    ip = inputParser;
    ip.addParameter('whittaker_lambda',1e7); 
    ip.addParameter('whittaker_d',2); 
    ip.addParameter('z_thresh',3); 
    ip.addParameter('spike_buffer',4); 
    ip.parse(varargin{:});
    for j=fields(ip.Results)'
        eval([j{1} '=ip.Results.' j{1} ';']);
    end

    % for a neutral reference: Whittaker–Henderson smoothing
    F_b = cell2mat(arrayfun(@(x) ...
        whittaker_smooth(F(:,x),'lambda',whittaker_lambda,'d',whittaker_d),...
        1:size(F,2), 'UniformOutput', false));
    
    % residual, then robust z-score (median/MAD, vs mean/SD)
    resid = F - F_b;
    medR = median(resid);
    madR = median(abs(resid - medR)) * 1.4826; % makes MAD match normal distrib st dev
    if madR == 0
        madR = eps;
    end
    robust_z = (resid - medR) ./ madR;
    
    % spikes to remove, dilated to include some data on either size
    spikes_to_rm = abs(robust_z) > z_thresh;
    spikes_to_rm = movmax(double(spikes_to_rm),[spike_buffer spike_buffer])>0;
        
    % nan 
    F_despiked = F;
    F_despiked(spikes_to_rm) = nan;    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% whittaker_smooth
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = whittaker_smooth(x, varargin)
    % y: raw signal, column vector
    % lambda: smoothing parameter (larger = smoother)
    % d: order of differences (typically 2 for second derivative)
    
    %%%  parse optional inputs %%%
    ip = inputParser;
    ip.addParameter('lambda',1e7); 
    ip.addParameter('d',2); 
    ip.parse(varargin{:});
    for j=fields(ip.Results)'
        eval([j{1} '=ip.Results.' j{1} ';']);
    end
    
    m = length(x);
    E = speye(m);
    D = diff(E, d);
    
    % Solve the penalized least squares system
    y = (E + lambda * (D' * D)) \ x;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FtoFc
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ Fc, scale, center ] = FtoFc( F, varargin )
    %FTOFC Normalizes by the 8th percentile in a sliding window and subtracts
    %   the median to create the modified DFoF
    
    %%%  parse optional inputs %%%
    ip = inputParser;
    ip.addParameter('scale_window',270); 
    ip.addParameter('baseline_prctile',0.08); 
    ip.addParameter('F_to_baseline',F); 
    ip.addParameter('spike_buffer',4); 
    ip.parse(varargin{:});
    for j=fields(ip.Results)'
        eval([j{1} '=ip.Results.' j{1} ';']);
    end
    
    
    scale = zeros(size(F));
    for i=1:size(scale,1)
        scale(i,:) = quantile(F_to_baseline(max(i-scale_window,1):min(i+scale_window,size(F,1)),:),baseline_prctile);
    end
    Fc = F./scale;
    center = median(Fc);
    Fc = bsxfun(@minus,Fc,center);

end




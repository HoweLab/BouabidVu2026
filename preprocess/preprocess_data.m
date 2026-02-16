%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preprocess_data
%%%%%%%%%%%%%%%%%%%%%%%%%%%

function output = preprocess_data(path_roi1,path_roi2,path_behav1,path_behav2,varargin)
    % takes as input paths to corresponding roi files and behav files from
    % a single recording. if only a 1-channel recording, leave path_roi2
    % and path_behav2 empty ([]).
    
    %%%  parse optional inputs %%%
    ip = inputParser;
    ip.addParameter('rm_neg_artifact',[0 0]); % negative going artifacts to remove
    ip.addParameter('rm_pos_artifact',[0 0]); % positive going artifacts to remove
    ip.addParameter('hp_hz',[0.3 0.1]); % high pass filter freq: default ch1 = 470nm = 0.3Hz, ch2 = 570nm = 0.1Hz
    ip.parse(varargin{:});
    for j=fields(ip.Results)'
        eval([j{1} '=ip.Results.' j{1} ';']);
    end
    
    
    % 1. align 2 channels if necessary
    if isempty(path_roi2) || isempty(path_behav2)
        output = struct;
        output.roi1 = load(path_roi1);
        output.behav1 = load(path_behav1);
    else
        output = align_2ch(path_roi1,path_roi2,path_behav1,path_behav2);        
    end        
    
    % roi deltaF/F preprocessing    
    for w = 1:sum(startsWith(fieldnames(output),'behav'))
        roiF = output.(['roi' num2str(w)]).F;
        sr = round(1/nanmean(diff(output.(['behav' num2str(w)]).timestamp)));
        
%         % sliding 8th percentile baseline
%         [Fc,F_baseline,~] = FtoFc(roiF,540);            
%         output.(['roi' num2str(w)]).FtoFcWindow = sr*30;
%         output.(['roi' num2str(w)]).Fc = Fc;
%         output.(['roi' num2str(w)]).F_baseline = F_baseline;

        % 2. calculate DFF from exponential baseline                                   
        [Fc,F_baseline] = FtoFc_exp(roiF,'exp_model','exp2');
        output.(['roi' num2str(w)]).F_baseline_exp = F_baseline;  
        output.(['roi' num2str(w)]).Fc_exp = Fc;

        % 3. highpass filter
        hp_cutoff = hp_hz(w);
        hp_steepness = 0.8;
        try
            output.(['roi' num2str(w)]).Fc_exp_hp = highpass(output.(['roi' num2str(w)]).Fc_exp,hp_cutoff,sr,'ImpulseResponse','fir','steepness',hp_steepness);
        catch exception
            disp('     could not high-pass filter')
            output.(['roi' num2str(w)]).Fc_exp_hp = [];
        end  
        output.(['roi' num2str(w)]).hp_steepness = hp_steepness;
        output.(['roi' num2str(w)]).hp_hz = hp_hz(w);

        % 4. if necessary, artifact removal 
        if rm_neg_artifact(w) == 1 || rm_pos_artifact(w) == 1
            try
                % use non-hp Fc
                tmp = remove_artifact(output.(['roi' num2str(w)]).Fc_exp,sr);
                tmp_roi = tmp.roi;
                % if there are negative artifacts to nan out
                if rm_neg_artifact(w) == 1
                    tmp_roi(tmp.neg_artifacts==1) = nan;
                end
                % if there are positive artifacts to nan out
                if rm_pos_artifact(w) == 1
                    tmp_roi(tmp.pos_artifacts==1) = nan;
                end                         
            catch exception
                disp('     could not remove artifact')
                tmp_roi = [];
            end
            % replace with nans
            if ~isempty(output.(['roi' num2str(w)]).Fc_exp_hp)                
                if ~isempty(tmp_roi)
                    output.(['roi' num2str(w)]).Fc_exp_hp_art = output.(['roi' num2str(w)]).Fc_exp_hp;
                    output.(['roi' num2str(w)]).Fc_exp_hp_art(isnan(tmp_roi)) = nan;
                else
                    output.(['roi' num2str(w)]).Fc_exp_hp_art = [];
                end
            end                      
        end  
        
        % 5. whether or not that day has signal, based on # transients
        % detected > chance
        if isfield(output.(['roi' num2str(w)]),'Fc_exp_hp_art') 
            if ~isempty(output.(['roi' num2str(w)]).Fc_exp_hp_art)
                output.(['roi' num2str(w)]).sig = get_signal_tr(output.(['roi' num2str(w)]).Fc_exp_hp_art);
            else
                output.(['roi' num2str(w)]).sig = 0;
            end
        elseif ~isempty(output.(['roi' num2str(w)]).Fc_exp_hp)
            output.(['roi' num2str(w)]).sig = get_signal_tr(output.(['roi' num2str(w)]).Fc_exp_hp);
        else
            output.(['roi' num2str(w)]).sig = 0;
        end        
    end
end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% align_2ch
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output = align_2ch(path_roi1,path_roi2,path_behav1,path_behav2)
    % Take as input paths of roi data files corresponding to ttlIn1 and 
    % ttlIn2 NIDAQ inputs to behavior file (e.g., 470nm, 570nm), behav 
    % files already aligned to TTLs (ttl1, ttl2). Aligns them based on 
    % timestamp, truncate and ignore non-simultaneous frames.
    %
    % Output a struct
    %
    % modified from codeo by Liangzhu Zhang

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
    output.roi1 = roi1;
    output.roi2 = roi2;
    output.behav1 = behav1;
    output.behav2 = behav2;    
    output.roi1idx = [I1start I1end];
    output.roi2idx = [I2start I2end];
    
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
% remove_artifact
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output = remove_artifact(roi,sr,varargin)
    % function roi = remove_artifact(roi)
    %
    % use a combination filter or DF/F (Fc) amplitude, .5Hz low-pass-filtered 
    % DF/F amplitude, and a 2s moving mean amplitude
    %
    % use a stringent filter to find periods to remove (-2.5*std)
    %
    % use a less-stringent filter to find the endpoints of those periods
    % (-2*std)
    %
    % Mai-Anh
    % 6/6/24


    %%%  parse optional inputs %%%
    ip = inputParser;
    ip.addParameter('neg_f1',-2);
    ip.addParameter('neg_f2',-2.5);
    ip.addParameter('pos_f1',2);
    ip.addParameter('pos_f2',5); % usually positive artifacts are huge
    ip.parse(varargin{:});
    for j=fields(ip.Results)'
        eval([j{1} '=ip.Results.' j{1} ';']);
    end

    % clunky but let's just do this to be simple: loop over ROIs
    output.roi = roi;
    output.roi_filtered = nan(size(roi));
    output.neg_artifacts = nan(size(roi));
    output.pos_artifacts = nan(size(roi));

    for r = 1:size(roi,2)
        % DF/F (Fc)
        fc = roi(:,r);
        % low-pass filtered
        f1p = lowpass(fc,.5,sr);
        % moving mean
        m = movmean(f1p,sr/2);

        % neg artifact: filters
        filter_1 = fc < neg_f1*std(fc) | f1p < neg_f1*std(f1p) | m < neg_f1*std(m);
        filter_2 = fc < neg_f2*std(fc) | f1p < neg_f2*std(f1p) | m < neg_f2*std(m);
        filter_3 = double(movmean(filter_2,sr*2)>0).*filter_1;
        output.neg_artifacts(:,r) = filter_3;

        % pos artifact: filters
        filter_1 = fc > pos_f1*std(fc) | f1p > pos_f1*std(f1p) | m > pos_f1*std(m);
        filter_2 = fc > pos_f2*std(fc) | f1p > pos_f2*std(f1p) | m > pos_f2*std(m);
        filter_3 = double(movmean(filter_2,sr*2)>0).*filter_1;
        output.pos_artifacts(:,r) = filter_3;

        % filtered ROI output
        fc_filtered = fc;
        fc_filtered(output.pos_artifacts(:,r)==1 | output.neg_artifacts(:,r)==1) = nan;    
        output.roi_filtered(:,r) = fc_filtered;    
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


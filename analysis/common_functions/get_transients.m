% function output = get_transients(roi)
%
% this function takes calculates a null distribution against which to
% compare to identify transients
%
% This combines a magnitude threshold with a duration threshold (e.g., 3
% timepoints must surpass a certain amplitude). The idea is bolster the
% amplitude threshold by ignoring spuriously high or low values, which are 
% unlikely to stay elevated (or depressed) for 3 (or however many, as 
% specified by 'n_surpass') frames.
%
% The bootstrapping of the null distribution works as follows: let's say, 
% you have an ROI timeseries that is 5000 elements long. This will randomly
% generate 10k (or n_bootstrap) timeseries of 5000 elements long by 
% randomly selecting from the input timeseries. This preserves the overall 
% magnitude and variance of the data but disrupts the temporal order, and,
% consequently, any auto-correlation. Then for each of the 10000
% timeseries, it calculates a sliding 3-pt (or # pts specified by 
% 'n_surpass') mean, and from these 10k (or n_it) resulting traces,
% calculates the percentiles 95:.5:100 (i.e., 95, 95.5, 96, ... 99.5, 100).
%
% Note that this takes the absolute values of magnitude, not making any
% assumptions about sign. So it tests whether a transient is significantly
% above or below the magnitude threshold as determined by the specified
% alpha value and generated null distribution.
%
% required input(s):
% roi               the m x n MATRIX of roi DFF timeseries where 
%                   m = timepoints and n = ROIs 
%                   OR the output STRUCT from this function.
%                   (Don't feed it the roi struct.)
%
% optional input(s):
% n_bootstrap       the number of bootstrapped samples (default 10k)
% n_surpass         the number of timepoints/frames/samples required to
%                   surpass the threshold to be considered a significant
%                   transient (default = 3)
% percentiles_lower the lower percentiles to calculate from the null distribution
%                   default = 0:.5:5 (0-5 in .5 increments)
% percentiles_upper the upper percentiles to calculate from the null distribution
%                   default = 95:.5:100 (95-100 in .5 increments)
% sig_alpha         the alpha for calculating significance. Note, this 
%                   assumes a one-tailed hypothesis test. For example, if
%                   you set alpha to 0.05, this calculates the upper
%                   significance using the 95th percentile and lower using
%                   the 5th percentile. For a two-tailed test, use alpha/2,
%                   e.g., alpha = 0.025, for an effective two-tailed test 
%                   with effective alpha of 5%. (default 0.025)
%
% output: the output struct
% This will contain all of the inputs (and the specified or default optional inputs)
% roi               the input roi DFF matrix
% n_bootstrap       the # of iterations for the bootstrap null
% n_surpass         the # of consecutive timepoints required to exceed null
% percentiles_lower list of lower bound percentiles calculated
% percentiles_upper list of upper bound percentiles calculated
% sig_alpha         the sig_alpha used (see above)
% null              the bootstrapped null distribution of the 3pt average
% null_upper        the upper percentile values of the null distribution
% null_lower        the lower percentiles values of the null distribution
% roi_sig_tr        same size as roi, 0 where no transients, numbered where
%                   there are transients. For example, 1st pos transient 
%                   gets 1s for the duration of the transient, 2nd gets 2s, 
%                   etc. Same for neg transients but (-) sign.
% pos_transients    struct with timestamps for onsets, offsets, peaks of 
%                   (+) transients, and magnitude of peaks
% neg_transients    struct with timestamps for onsets, offsets, peaks of 
%                   (-) transients, and magnitude of peaks

%                   
% Note that if you change your mind about the alpha value, you can easily
% recalculate roi_sig_pos by re-running this on the output struct and
% feeding in new arguments for sig_alpha.
%
% Mai-Anh Vu, updated 6/14/24


function output = get_transients(roi,varargin)

%%%  parse optional inputs %%%
ip = inputParser;
ip.addParameter('n_bootstrap',10000);
ip.addParameter('n_surpass',3);
ip.addParameter('percentiles_lower',0:0.5:5);
ip.addParameter('percentiles_upper',95:0.5:100);
ip.addParameter('sig_alpha',0.025);
ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end

% if the input was the roi DFF matrix, we need to calculate the null
% distribution percentiles

if isnumeric(roi) 
    output = struct;    
    output.roi = roi;       
    
    % calculate null distribution: get 10k samples of size n_surpass and
    % the mean of those 3    
    non_nan_idx = find(sum(isnan(roi),2) == 0);
    rand_idx = reshape(non_nan_idx(...
        randsample(numel(non_nan_idx),n_bootstrap*n_surpass,true)),...
        n_bootstrap,1,n_surpass);
    rand_idx_all = repmat(rand_idx,1,size(roi,2));
    roi_idx_all = repmat(1:size(roi,2),n_bootstrap,1,n_surpass);
    roi_null = reshape(roi(sub2ind(size(roi),rand_idx_all(:),roi_idx_all(:))),...
        size(rand_idx_all,1),size(rand_idx_all,2),size(rand_idx_all,3));
    
    % fill in some parameters
    output.n_bootstrap = n_bootstrap;
    output.n_surpass = n_surpass;
    output.percentiles_lower = percentiles_lower;
    output.percentiles_upper = percentiles_upper;
    output.sig_alpha = sig_alpha;    
    
    % now null distribution
    output.null = nanmean(roi_null,3); % this gets a null distribution for a 3-point event
    output.null_lower = prctile(output.null,percentiles_lower);
    output.null_upper = prctile(output.null,percentiles_upper);    
else
    output = roi;
end


% now calculate significance
prctile_idx_upper = find(output.percentiles_upper==round(200*(1-output.sig_alpha))/2); % find closest prctile (nearest .5) one-tailed-alpha
prctile_idx_lower = find(output.percentiles_lower==round(200*output.sig_alpha)/2); % find closest prctile (nearest .5) one-tailed-alpha
roi_sig_pos = double(output.roi > output.null_upper(prctile_idx_upper,:));
roi_sig_neg = double(output.roi < output.null_lower(prctile_idx_lower,:));
% account for nans (and don't count them as 0s)
roi_sig_pos(isnan(output.roi)) = nan;
roi_sig_neg(isnan(output.roi)) = nan;
% since we require n_surpass points to pass, ignore all values that don't
% surpass for n_surpass consecutive points, i.e., get rid of significant 
% points surpassing threshold if it's for less than n_surpass timepoints
for r = 1:size(output.roi,2)
    for i = 1:(n_surpass-1) 
        tmp = strfind(roi_sig_pos(:,r)',[0 ones(1,i) 0]);
        if ~isempty(tmp)
            tmp = repmat(tmp,(i+2),1) + repmat(transpose(0:(i+1)),1,numel(tmp));
            roi_sig_pos(tmp(:),r) = 0;                
        end
        tmp = strfind(roi_sig_neg(:,r)',[0 ones(1,i) 0]);
        if ~isempty(tmp)
            tmp = repmat(tmp,(i+2),1) + repmat(transpose(0:(i+1)),1,numel(tmp));
            roi_sig_neg(tmp(:),r) = 0;
        end
    end
end
roi_sig = roi_sig_pos + -1*roi_sig_neg;
output.roi_sig_tr = zeros(size(roi_sig)); % this puts the index of the significant transients where it exists

% now onsets and offsets of transients
signs = {'pos','neg'};
for s = 1:numel(signs)
    output.([signs{s} '_transients']).onset = cell(size(output.roi,2),1);
    output.([signs{s} '_transients']).offset = cell(size(output.roi,2),1);
    output.([signs{s} '_transients']).peak = cell(size(output.roi,2),1);
    output.([signs{s} '_transients']).magnitude = cell(size(output.roi,2),1);
end

signs_m = [1 -1];
for s = 1:numel(signs)
    for r = 1:size(output.roi,2)
        onsets = strfind(roi_sig(:,r)'==signs_m(s),[0 1]) + 1;
        offsets = strfind(roi_sig(:,r)'==signs_m(s),[1 0]);   
        if ~isempty(onsets) && ~isempty(offsets)
            if offsets(1) < onsets(1)
                offsets = offsets(2:end);
            end
            if onsets(end) > offsets(end)
                onsets = onsets(1:end-1);
            end            
        end
        output.([signs{s} '_transients']).onset{r} = vec(onsets);
        output.([signs{s} '_transients']).offset{r} = vec(offsets);
    end
end

% now peaks and troughs
for r = 1:size(output.roi,2)
    
    % peaks
    peaks = nan(size(output.pos_transients.onset{r},1),2);
    for i = 1:numel(output.pos_transients.onset{r})
        this_transient = output.roi(output.pos_transients.onset{r}(i):output.pos_transients.offset{r}(i),r);
        [amp,idx] = nanmax(this_transient);
        peaks(i,:) = [output.pos_transients.onset{r}(i)+idx-1 amp];
        output.roi_sig_tr(output.pos_transients.onset{r}(i):output.pos_transients.offset{r}(i),r) = i;
    end
    output.pos_transients.peak{r} = peaks(:,1);
    output.pos_transients.magnitude{r} = peaks(:,2);
    % troughs
    troughs = nan(size(output.neg_transients.onset{r},1),2);
    for i = 1:numel(output.neg_transients.onset{r})
        this_transient = output.roi(output.neg_transients.onset{r}(i):output.neg_transients.offset{r}(i),r);
        [amp,idx] = nanmin(this_transient);
        troughs(i,:) = [output.neg_transients.onset{r}(i)+idx-1 amp];
        output.roi_sig_tr(output.neg_transients.onset{r}(i):output.neg_transients.offset{r}(i),r) = -i;
    end
    output.neg_transients.peak{r} = troughs(:,1);
    output.neg_transients.magnitude{r} = troughs(:,2);
end
end
        
        
    
    


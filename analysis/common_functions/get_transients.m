% function output = get_transients(roi)
% Mai-Anh Vu, 9/16/2026


function output = get_transients(roi,pos_mad_k,neg_mad_k,varargin)

%%%  parse optional inputs %%%
ip = inputParser;
ip.addParameter('mad_clip_thresh',3);   % clipping bound, in MAD units (separate from detection k)
ip.addParameter('mad_max_iter',10);     % cap on clipping iterations
ip.addParameter('mad_tol',0.01);        % relative MAD change to declare convergence
ip.addParameter('tr_n',3);              % #timepoints in
ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end

% output
output = struct;
output.roi = roi;

% now do a sigma clipping procedure to estimate the noise MAD
[this_med,this_mad,this_frac,this_iter] = sigma_clip_mad(roi,...
    'clip_thresh',mad_clip_thresh,'max_iter',mad_max_iter,'tol',mad_tol);
output.mad.med = this_med;
output.mad.mad = this_mad;
output.mad.frac_clipped = this_frac;
output.mad.n_iter = this_iter;
    
% now find when the roi exceeds multipler x mad, and preserve nan
roi_sig_pos = roi > pos_mad_k * this_mad;
roi_sig_neg = roi < -neg_mad_k * this_mad;
% keep only exceedances meeting duration threshold
for r = 1:size(roi,2)
    for n = 1:(tr_n-1)        
        tmp = strfind(roi_sig_pos(:,r)',[0 ones(1,n) 0]);
        if ~isempty(tmp)
            tmp = repmat(tmp,(n+2),1) + repmat(transpose(0:(n+1)),1,numel(tmp));
            roi_sig_pos(tmp(:),r) = 0;                
        end
        tmp = strfind(roi_sig_neg(:,r)',[0 ones(1,n) 0]);
        if ~isempty(tmp)
            tmp = repmat(tmp,(n+2),1) + repmat(transpose(0:(n+1)),1,numel(tmp));
            roi_sig_neg(tmp(:),r) = 0;
        end
    end
end
roi_sig = roi_sig_pos + -1*roi_sig_neg;

% this is where i'll the index (and sign) sof each of the the significant transients where it exists
output.transients.sig = zeros(size(roi_sig)); 
% onsets and offsets of transients
signs = {'pos','neg'};
signs_m = [1 -1];
for s = 1:numel(signs)
    output.transients.(signs{s}).onset = cell(size(output.roi,2),1);
    output.transients.(signs{s}).offset = cell(size(output.roi,2),1);
    output.transients.(signs{s}).peak = cell(size(output.roi,2),1);
    output.transients.(signs{s}).magnitude = cell(size(output.roi,2),1);
end
for s = 1:numel(signs)
    for r = 1:size(output.roi,2)
        onsets = strfind(roi_sig(:,r)'==signs_m(s),[0 1]) + 1;
        offsets = strfind(roi_sig(:,r)'==signs_m(s),[1 0]);   
        if ~isempty(onsets) && ~isempty(offsets)
            offsets = offsets(offsets>onsets(1));
            if ~isempty(offsets)
                onsets = onsets(onsets<offsets(end));
            end
        end
        if ~isempty(onsets) && ~isempty(offsets)
            output.transients.(signs{s}).onset{r} = vec(onsets);
            output.transients.(signs{s}).offset{r} = vec(offsets);
        end
    end
end

% peaks and troughs
for r = 1:size(output.roi,2)
    % peaks
    peaks = nan(size(output.transients.pos.onset{r},1),2);
    for i = 1:numel(output.transients.pos.onset{r})
        this_transient = output.roi(output.transients.pos.onset{r}(i):output.transients.pos.offset{r}(i),r);
        [amp,idx] = nanmax(this_transient);
        peaks(i,:) = [output.transients.pos.onset{r}(i)+idx-1 amp];
        output.transients.sig(output.transients.pos.onset{r}(i):output.transients.pos.offset{r}(i),r) = i;
    end
    output.transients.pos.peak{r} = peaks(:,1);
    output.transients.pos.magnitude{r} = peaks(:,2);
    % troughs
    troughs = nan(size(output.transients.neg.onset{r},1),2);
    for i = 1:numel(output.transients.neg.onset{r})
        this_transient = output.roi(output.transients.neg.onset{r}(i):output.transients.neg.offset{r}(i),r);
        [amp,idx] = nanmin(this_transient);
        troughs(i,:) = [output.transients.neg.onset{r}(i)+idx-1 amp];
        output.transients.sig(output.transients.neg.onset{r}(i):output.transients.neg.offset{r}(i),r) = -i;
    end
    output.transients.neg.peak{r} = troughs(:,1);
    output.transients.neg.magnitude{r} = troughs(:,2);
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sensor_sigma_clip_mad: iteratively estimate a contamination-robust
% median/MAD for real-sensor fibers, excluding candidate real transients
% so they don't inflate the noise-scale estimate
%
% sigma_clip_mad: iteratively exclude candidate outliers (|x-med|>c*MAD)
% and recompute median/MAD on the remainder, per column (ROI), until
% convergence or max_iter is reached
function [med_clean,mad_clean,frac_clipped,n_iter] = sigma_clip_mad(x,varargin)

%%%  parse optional inputs %%%
ip = inputParser;
ip.addParameter('clip_thresh',3);
ip.addParameter('max_iter',10);
ip.addParameter('tol',0.01);
ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end

n_roi = size(x,2);
med_clean = nan(1,n_roi);
mad_clean = nan(1,n_roi);
frac_clipped = nan(1,n_roi);
n_iter = nan(1,n_roi);

for r = 1:n_roi
    this_x = x(:,r);
    this_x = this_x(~isnan(this_x));  % drop pre-existing artifact-masked samples
    n_valid = numel(this_x);
    included = true(n_valid,1);
    prev_mad = nan;
    this_med = median(this_x);
    this_mad = mad(this_x,1);
    for it = 1:max_iter
        this_med = median(this_x(included));
        this_mad = mad(this_x(included),1);
        if ~isnan(prev_mad) && this_mad>0 && abs(this_mad-prev_mad)/prev_mad < tol
            break
        end
        prev_mad = this_mad;
        included = abs(this_x-this_med) <= clip_thresh*this_mad;
    end
    med_clean(r) = this_med;
    mad_clean(r) = this_mad;
    frac_clipped(r) = 1-sum(included)/n_valid;
    n_iter(r) = it;
end
end
    


% function output = local_morans_I_sig_moran(moran_map,null_moran_struct,how_best,vol_thresh)
%
% returns the hottest hotspot from comparison of local moran's with null moran
%
% inputs:
% moran_map: map of local_morans_I
% null_moran: output from local_morans_I_bootstrap_null.m
% how_best = how to determine the hottest hotspot: 'dominant','min',or,'max'
%
% optional inputs:
% vol_thresh: minimum volume for hotspot
% generate_rand: the number of random connected volumes to generate of same 
%       size as hotspot (set to 0 to not do this)
% mask: optional mask to mask moran map and also from within which to
%       generate random volumes (if opt in)
%
% Mai-Anh Vu, 2025
function output = local_morans_I_sig_moran(smooth_map,moran_map,null_moran,how_best,varargin)

%%%  parse optional inputs %%%
ip = inputParser;
ip.addParameter('vol_thresh',10^3);
ip.addParameter('generate_rand',10000); % set to 0 to not do this
ip.addParameter('mask',[])
ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end
if isempty(mask)
    mask = ones(size(moran_map));
end

% significant moran
sig_map = moran_map > null_moran.prctile95 &  mask==1;

% connected components
cc = bwconncomp(sig_map);
cc_size = cellfun(@(x) numel(x),cc.PixelIdxList);
cc_vals = cellfun(@(x) mean(smooth_map(x)),cc.PixelIdxList);
keep_idx = cc_size>vol_thresh; % more than 500um^3
cc_size = cc_size(keep_idx);
cc_vals = cc_vals(keep_idx);
cc.PixelIdxList = cc.PixelIdxList(keep_idx);
% sorting, depending on condition
if strcmp(how_best,'max')
    [~,sort_idx] = sort(cc_vals,'descend');
elseif strcmp(how_best,'min')
    [~,sort_idx] = sort(cc_vals,'ascend');
elseif strcmp(how_best,'dominant')
    [~,sort_idx] = sort(abs(cc_vals),'descend');
end

% make a mask of this volume
sig_map2 = zeros(size(sig_map));
sig_map2(cc.PixelIdxList{sort_idx(1)}) = 1; 
output.map = sig_map2;
output.vox = cc.PixelIdxList{sort_idx(1)};

% also generate a random distribution of volumes of this size
if generate_rand > 0
    output.rand =  rand_volume(size(mask),numel(output.vox),'n',generate_rand,'mask',mask);
end

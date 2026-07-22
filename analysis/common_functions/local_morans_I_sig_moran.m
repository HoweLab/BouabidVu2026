% function output = local_morans_I_sig_moran(moran_struct,null_moran_path,varargin)
%
% returns the hottest hotspot from comparison of local moran's with null moran
%
% inputs:
% moran_map: map of local_morans_I
% null_moran: output from local_morans_I_bootstrap_null.m save_null option
%
% optional inputs:
% generate_rand: the number of random connected volumes to generate of same 
%       size as hotspot (set to 0 to not do this)
% mask: optional mask to mask moran map and also from within which to
%       generate random volumes (if opt in)
%
% Mai-Anh Vu, 2025
% edited 2026/07/22 to incorporate boostrapped null cluster size

function moran_struct = local_morans_I_sig_moran(moran_struct,null_moran_path,varargin)

%%%  parse optional inputs %%%
ip = inputParser;
ip.addParameter('mask',[])
ip.addParameter('alpha_val',0.05)
ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end

% mask
if isempty(mask)
    mask = ones(size(moran_map));
end

% corresponding percentile
this_prctile = 100*(1-alpha_val);
p_field = (['prctile' strrep(sprintf('%2.02f',this_prctile),'.','p')]);
vox_thresh = moran_struct.null_stats.(p_field);


% maximum cluster at each iteration of null moran
null_moran = matfile(null_moran_path); 
n_it = size(null_moran.null_moran,4);
null_max_clusters = nan(n_it,1);

for i = 1:n_it
    this_null = null_moran.null_moran(:,:,:,i);
    this_null(mask==0) = nan;
    this_null_sig = this_null > vox_thresh;
    sig_clusters = bwconncomp(this_null_sig);
    sig_clusters = cellfun(@numel, sig_clusters.PixelIdxList);
    null_max_clusters(i) = max(sig_clusters);
end
% save out cluster threshold
clust_thresh = prctile(null_max_clusters,this_prctile);
moran_struct.null_stats.(strrep(p_field,'prctile','clustsize')) = clust_thresh;
    
% actual significance map
sig_map = moran_struct.moran;
sig_map(mask==0) = nan;
sig_map = sig_map > vox_thresh;

% hotspots
sig_clusters = bwconncomp(sig_map);
sig_clusters_size = cellfun(@(x) numel(x),sig_clusters.PixelIdxList);
keep_idx = sig_clusters_size > clust_thresh; 
moran_struct.sig.hotspot = sig_clusters.PixelIdxList(keep_idx);

% also generate a random distribution of volumes of hotspot(s)
if generate_rand > 0
    for i = 1:numel(moran_struct.sig.hotspot)
        moran_struct.sig.rand{i} = rand_volume(size(mask),...
            numel(moran_struct.sig.hotspot{i}),...
            'n',generate_rand,'mask',mask);
    end
end

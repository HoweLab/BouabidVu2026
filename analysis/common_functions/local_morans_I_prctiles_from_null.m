% function moran_struct = local_morans_I_prctiles_from_null(moran_struct,null_moran_path,prctiles)
%
% While calculation of voxelwise prctiles is already built into
% local_morans_I_boostrap_null, this function enables it post-hoc
% 
% inputs:
% moran_struct      - the output struct from local_moran_I_bootstrap_null
% null_moran_path   - output from local_morans_I_bootstrap_null save_null option
% prctiles          - array of percentiles to calculate (e.g., [95 97.5 99 99.5 99.9])
%
% returns:
% moran_struct with updated prctiles field



function moran_struct = local_morans_I_prctiles_from_null(moran_struct,null_moran_path,prctiles)

% the null moran mat
null_moran = matfile(null_moran_path);
n_it = size(null_moran.null_moran,4);

% keep an buffer of upper values to be able to calculate prctiles
n_to_keep =  ceil((100-min(prctiles))/100*n_it);
upper_buffer = [];

% read in an iteration at at time, keeping only as many as necessary to
% calculate the desired prctiles
for i = 1:n_it
    this_null = null_moran.null_moran(:,:,:,i);
    upper_buffer = cat(4,upper_buffer,this_null);
    upper_buffer = sort(upper_buffer,4,'descend');
    upper_buffer = upper_buffer(:,:,:,1:n_to_keep);
end
    
% output prctiles
upper_buffer = sort(upper_buffer,4,'ascend');
if ~isfield(moran_struct,'null_stats')
    moran_struct.null_stats = struct;
end
for p = 1:numel(prctiles)
    p_field = (['prctile' strrep(sprintf('%2.02f',prctiles(p)),'.','p')]);
    p_idx = round(prctiles(p)/100*n_it-n_it+n_to_keep+1);
    moran_struct.null_stats.(p_field) = upper_buffer(:,:,:,p_idx);
end





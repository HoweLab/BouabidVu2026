% output = rand_volume(input_size,volume_size,varargin)
%
% randomly select a contiguous volume of specified size (volume_size) from
% within a larger volume of specieid size (input_size)
%
% inputs:
% input_size: 1 x 3 vector
% volume_size: a single integer = the number of voxels in the selected
%   volume
%
% optional
% mask: a mask from within which to select the volume - can be a volume or
%   a list of voxels
% n: the number of volumes desired (default = 1)
%
% output: a cell array of voxel lists
%
% Mai-Anh Vu, 9/17/2025, as a way to generate a null distribution to
% quantify overlap between "hotspot" volumes
%
%

function output = rand_volume(input_size,volume_size,varargin)

    %%%  parse optional inputs %%%
    ip = inputParser;
    ip.addParameter('mask',[]);
    ip.addParameter('n',1);
    ip.parse(varargin{:});
    for j=fields(ip.Results)'
        eval([j{1} '=ip.Results.' j{1} ';']);
    end

    % our input volume
    input_vol = zeros(input_size);
    if ~isempty(mask)
        input_vol(mask==1) = 1;
    else
        input_vol = input_vol + 1; % make it all ones
    end
    output = cell(n,1);
    parfor i = 1:n
        output_vol = get_rand_vol(input_vol,volume_size);
        output{i} = find(output_vol==1);
        if rem(i,100)==0 % display every 100
            disp(i)
        end
    end

end


% generate a contiguous connected volume of desired size from within an
% input mask
function output_vol = get_rand_vol(input_vol,volume_size)
    rng('shuffle')
    % output volume
    output_vol = zeros(size(input_vol));
    % input mask indices
    input_vox = find(input_vol==1);
    % start with a seed voxel
    output_vol(input_vox(randsample(numel(input_vox),1))) = 1;
    % Iterate: dilate by one voxel all around to generate candidate voxels
    % then select some random number of those voxels to add.
    % The number of voxels is a random integer between 1 and 1/2 the number
    % of candidate voxels. Do this until we've reached the desired volume.
    while sum(output_vol(:))<volume_size
        candidate_voxels = find((imdilate(output_vol,strel('sphere',1))-output_vol)==1);
        candidate_voxels = intersect(input_vox,candidate_voxels);
        n_samp = min([volume_size-sum(output_vol(:)) randi(floor(numel(candidate_voxels)/2))]);
        output_vol(candidate_voxels(randsample(numel(candidate_voxels),n_samp,'false'))) = 1;
    end
end
    

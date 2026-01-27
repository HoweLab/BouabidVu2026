% get indices for map volume from fiber table. 
% str, the output from get_striatum_vol_mask serves as the reference
%
% Mai-Anh Vu, 2025

function map_ind = get_map_ind(fib,str)
    map_idx =  nan(size(fib,1),3); % DV, ML, AP
    [~,map_idx(:,1)] = arrayfun(@(x) min(abs(x-str.info.DV)),fib.fiber_bottom_DV);
    [~,map_idx(:,2)] = arrayfun(@(x) min(abs(x-str.info.ML)),fib.fiber_bottom_ML);
    [~,map_idx(:,3)] = arrayfun(@(x) min(abs(x-str.info.AP)),fib.fiber_bottom_AP);
    map_ind = sub2ind(size(str.striatum_mask),map_idx(:,1),map_idx(:,2),map_idx(:,3));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_mouse_neighborhood_r:
%
% Based on a mouse's fiber locations, determine the minimum radius r such
% that for any fiber, there are x fibers within an r-radius sphere.
% 
% x can either be hard-coded, or determined via Kelejian and Prucha (2007),
% which says the cube-root of the number of fibers
%
% Mai-Anh Vu, 2026/07/15
%
function mouse_neighborhood = get_mouse_neighborhood_r(fib,varargin)
    %%%  parse optional inputs %%%
    ip = inputParser;
    ip.addParameter('min_n_fibs',[]);   % #fibers required to be within radius r    
    ip.parse(varargin{:});
    for j=fields(ip.Results)'
        eval([j{1} '=ip.Results.' j{1} ';']);
    end
    
    % initialize output
    mouse_neighborhood = struct;
    
    % make this a cell
    if ischar(fib.mouse)
        fib.mouse = cellstr(fib.mouse);
    end
    
    % loop over mice
    mice = unique(fib.mouse);
    for m = 1:numel(mice)
        mouse = mice{m};
        this_fib = fib(ismember(fib.mouse,(mouse)),:);
        if isempty(min_n_fibs)
            this_min_n_fibs = round(size(this_fib,1)^(1/3));
        else
            this_min_n_fibs = min_n_fibs;
        end
        
        % get fib-fib distances
        fib_dist = pdist2([this_fib.fiber_bottom_DV this_fib.fiber_bottom_ML this_fib.fiber_bottom_AP],...
            [this_fib.fiber_bottom_DV this_fib.fiber_bottom_ML this_fib.fiber_bottom_AP]);
        % now sort each column
        fib_dist = sort(fib_dist,1);
        dist_req = max(fib_dist(this_min_n_fibs,:));
        % output
        mouse_neighborhood.(mouse).r = dist_req;
        mouse_neighborhood.(mouse).n = this_min_n_fibs;
    end
end
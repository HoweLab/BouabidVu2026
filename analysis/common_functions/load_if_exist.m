% output = check_if_exist(input_path,varargin)
%
% checks if an output exists, and if so loads and returns it;
% otherwise, returns a blank struct or array

function output = load_if_exist(output_path,varargin)


%%%  parse optional inputs %%%
ip = inputParser;
ip.addParameter('output_class','struct'); % could also be other classes like 'double','cell',etc
ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end

if exist(output_path,'file')
    output = load(output_path);
else
    eval(['output = ' output_class '.empty']);
end
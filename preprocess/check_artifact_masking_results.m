% check_artifact_masking_results(mouse,expdate,varargin)
%
% plots results from classifyMultiChannelArtifacts as run via
% preprocess_data.m
%
% can also give it the data struct itself as the mouse argument; leave
% expdate empty
%
% optional inputs:
% data_dir      - expect data path to be data_dir\mouse\expdate\mouse_expdate.mat; default 'G:'
% channel_names - default {'ACh','DA'};
% roi_nums      - which ROIs to plot; default is [] which will plot all
% Fc_field      - which DFF field to plot; default is 'Fc_exp'
% n_rows        - each ROI will be plotted on one row; how many rows per figure
% artifact_mask - field name of the artifact mask
%
% Mai-Anh Vu, 2026/08/03
%

function check_artifact_masking_results(mouse,expdate,varargin)


%%%  parse optional inputs %%%
ip = inputParser;
ip.addParameter('data_dir','G:');                       % expect data_dir\mouse\expdate\mouse_expdate.mat
ip.addParameter('channel_names',{'ACh','DA'});          % channel_names; see preprocess_data.m
ip.addParameter('roi_nums',[]);                         % which ROIs to plot (default = all)
ip.addParameter('Fc_field','Fc');                       % which DFF field (default = 'Fc_exp')
ip.addParameter('art_mask_suffix','');                  % suffix of artifact mask fields (artifact1_suffix,artifact2_suffix, artifact3_suffix,artifact_mask_suffix)
ip.addParameter('art_colors',[1 0 1; 0 1 0; 0 1 1; 1 1 0]);    % default: type1 = magenta, type2 = green, type3 = cyan
ip.addParameter('art_markers',{'o','.','x'});            
ip.addParameter('art_marker_sizes',[3 6 3]);            
ip.addParameter('datafile_suffix','');                  % datafile format [mouse]_[expdir][datafile_suffix].mat; default '', assume format MOUSE_EXPDIR.mat 
ip.addParameter('n_rows',8);                          	% how many ROIs per figure
ip.addParameter('fig_pos',[]);                          % figure position; if blank, fullscreen
ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end

if isempty(expdate)
    data = mouse;
    mouse = [];
    expdate = [];
else
    if exist(fullfile(data_dir,mouse,expdate,[mouse '_' expdate datafile_suffix '.mat']),'file')
        data = load(fullfile(data_dir,mouse,expdate,[mouse '_' expdate datafile_suffix '.mat']));
    else
        disp('Data file not in expected location.')
        return
    end
end

if isempty(roi_nums)
    roi_nums = 1:size(data.(channel_names{1}).(Fc_field),2);
end
f = 0; % figure counter
k = 0; % row counter
for r = 1:numel(roi_nums)
    if k == 0
        f = f+1;
        if ~isempty(fig_pos)
            figure('Position',fig_pos)
        else
            figure('units','normalized','outerposition',[0 0 1 1],'visible','on');
        end
    end
    k = k + 1;
    this_row = k;
    for c = 1:numel(channel_names)
        this_data = data.(channel_names{c}).(Fc_field)(:,roi_nums(r));        
        this_col = c;
        subplot(n_rows,numel(channel_names),(this_row-1)*numel(channel_names)+this_col)
        hold on
        plot(this_data,'-k')
        for a = 1:3
            this_art = find(data.(channel_names{c}).(['artifact' num2str(a) art_mask_suffix])(:,roi_nums(r))==1);
            plot(this_art,this_data(this_art),art_markers{a},'Color',art_colors(a,:),'MarkerSize',3); 
        end
        title([mouse ': ' expdate ' | ' channel_names{c} ': ' num2str(roi_nums(r))])
        set(gca,'XLim',[1 size(this_data,1)])
    end
    
    if k == n_rows % reset if we're at the end
        k = 0;
    end
end
        
        
    

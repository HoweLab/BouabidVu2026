%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot_hotspot_sig_sites:
function plot_hotspot_sig_sites(this_map,str,sig,fib,varargin)

%%%  parse optional inputs %%%
ip = inputParser;
% directories and atlas: see https://github.com/HoweLab/MultifiberLocalization    
ip.addParameter('y_max',[]);
ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end

map_ind = get_map_ind(fib,str);
in_map = this_map(map_ind) == 1;

sig_in = sig(in_map==1);
sig_out = sig(in_map==0);

figure('Position',[100 100 450 300])
hold on

if ~isempty(sig_in)
    bar(1,numel(sig_in),'FaceColor',[.9 .9 .9],'HandleVisibility','off')
    bar(1,sum(sig_in),'FaceColor',[.5 .5 .5],'HandleVisibility','off')
end

if ~isempty(sig_out)
    bar(2,numel(sig_out),'FaceColor',[.9 .9 .9])
    bar(2,sum(sig_out),'FaceColor',[.5 .5 .5])
end
legend({'Insignificant','Significant'},'Location','NorthWest')
set(gca,'XTick',[1 2],'XTickLabel',{'In','Out'})
ylabel('# Sites')
if ~isempty(y_max)
    set(gca,'YLim',[0 y_max])
end
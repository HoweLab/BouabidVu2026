% Mai-Anh Vu, 2025
function plot_lag_corr_scatter(corr_r,corr_lat,varargin)

%%%  parse optional inputs %%%
ip = inputParser;
% directories and atlas: see https://github.com/HoweLab/MultifiberLocalization    
ip.addParameter('clust_id',[]);
ip.addParameter('r_bins',-1:0.1:1);
ip.addParameter('lat_bins',(-9/18*1000):(1000/18):(9/18*1000));
ip.addParameter('voxel_size_rounding',1); % e.g., if voxel size is 0.05mm, round such that voxel boundaries are multiples of 0.05
ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end
if isempty(clust_id)
    clust_id = ones(size(corr_r));
    cat_colors = [.5 .5 .5];
else
    cat_colors = lines(max(clust_id));
end


figure('Position',[300 300 500 450])

% corr coeff r
ax(1) = subplot(4,4,[1 2 3]);
hold on
for i = 1:max(clust_id)
    histogram(corr_r(clust_id==i),'BinEdges',r_bins,...
        'Normalization','probability','FaceColor',cat_colors(i,:))
end
xline(0)
set(gca,'XTickLabel',[],'FontSize',12)

% corr latencies
ax(2) = subplot(4,4,[8 12 16]);
hold on
for i = 1:max(clust_id)
    histogram(corr_lat(clust_id==i),'BinEdges',lat_bins,...
        'Normalization','probability','FaceColor',cat_colors(i,:),...
        'orientation','horizontal')
end
yline(0)
set(gca,'YTickLabel',[],'FontSize',12)

% scatter plot
ax(3) = subplot(4,4,[5 6 7 9 10 11 13 14 5]);
hold on
for i = 1:max(clust_id)
    plot(corr_r(clust_id==i),corr_lat(clust_id==i),'o',...
        'Color',cat_colors(i,:)*.5,'MarkerFaceColor',cat_colors(i,:))
end
xline(0)
yline(0)
linkaxes(ax([1 3]),'x')
linkaxes(ax([2 3]),'y')
set(gca,'YLim',[min(lat_bins) max(lat_bins)],...
    'XLim',[min(r_bins) max(r_bins)],'FontSize',12)
xlabel('Pearson \it{r}','interpreter','tex')
ylabel('Latency (ms)')
    

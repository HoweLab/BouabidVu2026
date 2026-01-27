% plot_violin_in_out
%
% Mai-Anh Vu, 2025
% str = the returned output from get_striatum_vol_mask.m 
% hotspot_map has to be the same size as str.striatum_mask

function plot_violin_in_out(vals,fib,str,hotspot_map)

% figure out which values are inside v outside hotspot map
map_ind = get_map_ind(fib,str);
in_map = hotspot_map(map_ind) == 1;
vals_in = vals(in_map==1);
vals_out = vals(in_map==0);

figure
hold on
violinPlotSingle(vals_in,'X',1,'violinColor',[0 0 0],'PointSize',10)
violinPlotSingle(vals_out,'X',2,'violinColor',[.5 .5 .5],'PointSize',10)
yline(0)
set(gca,'XTick',[1 2],'XTickLabel',{'In','Out'},'FontSize',12,'XLim',[0.5 2.5])
[p,~] = ranksum(vals_in,vals_out);
title([' p = ' sprintf('%0.3f',p)],'Interpreter','none')

end


function violinPlot(data,varargin)

% violinPlot(data,'OptionalParam1',val1,'OptionalParam2,val2,...)
% My own violin plot code based on https://github.com/bastibe/Violinplot-Matlab.
% Generates a figure with violin plots. Data should be n x p, where p is the number 
% of categories. Typical axis parameters can be edited with set(gca,'Param1',val1,...)
%
% 4/22/20 Mai-Anh Vu
%
%	 	
%	'Bandwidth'			Bandwidth of the kernel density estimate.
%						Should be between 10% and 50% of the data range
% 	'ViolinWidth'		Width of violin in axis space. Defaults to 0.3.
%	'ViolinColor'		n x 3 color matrix to define colors per violin. 
%						Default is Matlab "colors" colormap.
%	'ViolinAlpha'		Transparency of the violin area. Default 0.3.
%	'EdgeColor'			Color of violin outline. Default same as ViolinColor.
%	'EdgeWidth'			Width of violin outline. Default 2. Set 0 for no outline.
%	'PointSize'			Marker size for points. Default 10
%	'PointAlpha'		Transparency of the Points. Default 0.5
% 	'PointsOnTop'		Points on top (vs under) violin. Default 0;
%	'PointColor'		Color of points. Default ViolinColor.
%	'MedianColor'		Fill color of the median. Defaults to [1 1 1]
%	'MedianPointSize'	Marker size of the median. Defaults to 50;
%	'BoxColor'			Color of the whiskers and outline of the median point. 
%                       Default [.5 .5 .5]
% 	'BoxLineWidth'		Line width of the whiskers and outline of the median point. 
%                       Default 3.
% 	'BoxIQRLineWidth'	Line width of the extended whiskers (1.5 IQR). Default 1.
% 	'ShowData'			Whether to show data points. Default true (1)
%	'ShowMean'			Show mean. Default false (0).
%   'MeanOffset'        How far off from the midline to plot mean 
%                       (to avoid overlap with boxplot). Default 0.
%	'ShowSEM'			Show SEM (std/sqrt(n)). Default false (0).
%	'MeanColor'			Color of mean (and SEM). Default [0 0 0].
%	'MeanLineWidth'		Width of mean (and SEM) lines. Default 2.
%	'MeanWidth'			Width of mean line marker. Default 0.1.
%   'X'                 where to put the violin along the X-axis




%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  parse inputs %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
ip = inputParser;
ip.addParameter('Bandwidth',[]);
ip.addParameter('ViolinWidth',0.3);
ip.addParameter('ViolinColor',lines(1));
ip.addParameter('ViolinAlpha',0.3);
ip.addParameter('EdgeColor',[])
ip.addParameter('EdgeWidth',1);
ip.addParameter('ShowData',1);
ip.addParameter('PointSize',10);
ip.addParameter('PointAlpha',0.5);
ip.addParameter('PointsOnTop',0);
ip.addParameter('PointColor',[])
ip.addParameter('MedianColor',[1 1 1]);
ip.addParameter('MedianPointSize',75);
ip.addParameter('BoxColor',[.5 .5 .5]);
ip.addParameter('BoxLineWidth',3);
ip.addParameter('BoxIQRLineWidth',1);
ip.addParameter('ShowMean',0);
ip.addParameter('MeanOffset',0);
ip.addParameter('ShowSEM',0);
ip.addParameter('MeanColor',[0 0 0]);
ip.addParameter('MeanLineWidth',2);
ip.addParameter('MeanWidth',0.1);
ip.addParameter('x',1);

ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end
if isempty(EdgeColor)
    EdgeColor = ViolinColor;
end
if isempty(PointColor)
    PointColor = ViolinColor;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  main loop %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate kernel density estimation for the violin
[density, value] = ksdensity(data, 'bandwidth', Bandwidth);
bounds = [find(value<=min(data),1,'last') find(value>=max(data),1,'last')];
density = density(bounds(1):bounds(2));
value = value(bounds(1):bounds(2));	
% if all data is identical
if min(data) == max(data)
    density = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now onto plotting violin & points
jitterWidth = ViolinWidth/max(density);

% plot the points within the violin area if they're going under the violin
if ShowData==1 && PointsOnTop==0
    if numel(density)>1
        jitterstrength = interp1(value, density*jitterWidth, data);
    else
        jitterstrength = density*jitterWidth;
    end
    jitter = 2*(rand(size(data,1),1)-0.5);
    if size(PointColor,1)==1
        scatter(x+jitter.*jitterstrength,data,PointSize,...
            'MarkerFaceColor',PointColor,...
            'MarkerFaceAlpha',PointAlpha,...
            'MarkerEdgeColor','none');
    else
        for p = 1:numel(data)
            scatter(x+jitter(p).*jitterstrength(p),data(p),PointSize,...
            'MarkerFaceColor',PointColor(p,:),...
            'MarkerFaceAlpha',PointAlpha,...
            'MarkerEdgeColor','none');
        end
    end
end

% plot the violin
vx = [(x + density*jitterWidth) fliplr((x - density*jitterWidth))];
vy = [value fliplr(value)];
f = fill(vx,vy,ViolinColor,...
    'EdgeColor',EdgeColor,...		
    'FaceAlpha',ViolinAlpha);
if EdgeWidth>0
    set(f,'LineWidth',EdgeWidth);
else
    set(f,'LineStyle','none');
end

% if points are on top, plot the points now
if ShowData==1 && PointsOnTop==1
    if numel(density)>1
        jitterstrength = interp1(value, density*jitterWidth, data);
    else
        jitterstrength = density*jitterWidth;
    end
    jitter = 2*(rand(size(data,1),1)-0.5);
    if size(PointColor,1)==1
        scatter(x+jitter.*jitterstrength,data,PointSize,...
            'MarkerFaceColor',PointColor,...
            'MarkerFaceAlpha',PointAlpha,...
            'MarkerEdgeColor','none');
    else
        for p = 1:numel(data)
            scatter(x+jitter(p).*jitterstrength(p),data(p),PointSize,...
            'MarkerFaceColor',PointColor(p,:),...
            'MarkerFaceAlpha',PointAlpha,...
            'MarkerEdgeColor','none');
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now onto plotting some stats    

% plot 1.5IQR, quartiles,  and median
quartiles = quantile(data,[0.25 .5 0.75]);
IQR = quartiles(3) - quartiles(1);
extendedWhiskers = [quartiles(1)-1.5*IQR quartiles(3)+1.5*IQR];
plot([x x],extendedWhiskers,'-k',...
    'Color',BoxColor,...
    'LineWidth',BoxIQRLineWidth)
plot([x x],[quartiles(1) quartiles(3)],'-k',...
    'Color',BoxColor,...
    'LineWidth',BoxLineWidth);
scatter(x,quartiles(2),MedianPointSize,...
    'MarkerFaceColor',MedianColor,...
    'MarkerEdgeColor',BoxColor);


% plot the mean and sem, if applicable
if ShowMean==1
    meanValue = mean(data);
    plot([x-MeanOffset-MeanWidth x-MeanOffset+MeanWidth],[meanValue meanValue],'-k',...
        'Color',MeanColor,...
        'LineWidth',MeanLineWidth);
    if ShowSEM==1
        semValue = std(data)/sqrt(size(data,1));
        errorbar(x-MeanOffset,meanValue,semValue,...
            'Color',MeanColor,...
            'LineWidth',MeanLineWidth);
    end
end	

end
% plot_site_maps(values, fib, varargin) 
%
% This function makes a 3D heatmap circle plot of values from a grid 
% recording (x = ML, y = AP, z = DV). The first input is the data. The 
% second input can either be a calib.mat file or the output table from 
% fiber localization. The data has to be in the same order as the rows in 
% the calibration struct or fiber localization table. 
%
% The default view is a sagittal projection but you can change this manually 
% by swiveling the output, programmatically via Matlab's "view" command on 
% the output figure to set the azimuth and elevation angles, or 
% programmatically by using the 'viewAngle' optional input to this function
% to supply the azimuth and elevation angles. See Matlab's "view" command 
% for more information. Note that since this one is 3D, there's no scaling 
% based on depth as in bubblePlot. 
%
% There are a bunch of other possible optional inputs for customization. 
% Open the script to see more.
%
% Note You can also use a subset (e.g., you made a new table using a subset
% of the rows)  of the table from the fiber localization. Note, however,
% that the ROI field has to be updated to reflect the proper number of
% ROIs. For example:
% someRows = [1 2 3 5 6 8 10 30 40]; % the rows i want for my new table
% subTbl = fibLocTable(someRows,:); % take a subset of the table
% subTbl.ROI = transpose(1:size(subTbl,1)); % reset the ROI#
%
% Mai-Anh 8/12/222
% 
% updated 11/16/2022 to input colormap directly as 'colormapOption' 
%   - this can either be the string to a colormap (e.g., 'parula'), 
%   or a n x 3 matrix of colors
% updated 11/16/2022  # colorbar ticks can be set with 'colorbarTickDec'
% updated 5/17/2023 transparecy can be set for marker colors and outline
% updated 7/13/2023 to remove default axis lims 
% updated 7/28/2023 custom colormap fix
% updated 5/15/2024 to handle DV sign (whether we have it all positive
%   or negative)
%
function plot_site_maps(values, fib, varargin)
%%%  parse optional inputs %%%
ip = inputParser;
%%% basic settings
ip.addParameter('str',[]); % if striatum struct is supplied, will plot striatum outlines
ip.addParameter('saveDir',[]); % directory to save figure (if blank, won't save) 
%%% bubble appearance
ip.addParameter('cmapBounds',[min(values) max(values)]); % manually set the color bounds (note RedBlue will automatically make it symmetric)
ip.addParameter('outlineID',ones(size(values))); % the indices that go with outline color or width (i.e., if you have 2 colors, you can specify whether the fiber should be color 1 or 2)
ip.addParameter('outlineColor',[0 0 0]);% outline color: one row per color. row indices correspond with outlineID
ip.addParameter('outlineWidth',1);% outline width: you can have multiple widths corresponding to outlineID if you want
ip.addParameter('outlineWidthID',ones(size(values)));% outline width ID: which bubbles have which outlineWidth values
ip.addParameter('outlineAlpha',1);% outline alpha: you can have multiple transparencies corresponding to outlineAlphaID if you want
ip.addParameter('outlineAlphaID',ones(size(values)));% outline alpha ID: which bubbles have which outlineAlpha values
ip.addParameter('colorAlpha',1);% color alpha: you can have multiple bubble fill transparencies corresponding to colorAlphaID if you want
ip.addParameter('colorAlphaID',ones(size(values)));% color alpha ID: which bubbles have which colorAlpha values
ip.addParameter('blankBubble',[]); % the list of bubbles to leave blank
ip.addParameter('blankBubbleColor','none'); % the color for blank bubbles
ip.addParameter('useStrDepth',0); % use approx depth from striatal dorsal surface (instead of fiber length or DV) - NOTE: this option hasn't been fully tweaked. it's just approximate.
%%% bubble size options
ip.addParameter('bubbleSize',100); % size of bubbles
%%% color map and colorbar options
ip.addParameter('includeColorbar',1); % whether or not to include the colorbar
ip.addParameter('colorbarLoc','eastoutside'); % location of the colorbar (check out matlab 'colorbar' for more info)
ip.addParameter('colorbarTickDec',2); % number of tick marks on the colorbar
ip.addParameter('colormapOption','parula'); % colormap (either the string to a map, or a matrix)
ip.addParameter('colormapBins',256); % the number of color bins for the colormap
ip.addParameter('colorbarLabel',[]); % a text label for the colorbar
%%% axes options
ip.addParameter('newFig',0); % whether to make this in a new figure  (default = 1 = new fig, 0 is plot on current fig)
ip.addParameter('axesLabels',1); % whether to label axes (tick marks, numbers, etc): AP, ML, DV axes titles are labeled, regardless
ip.addParameter('projection','axial'); % default projection is axial; other options 'coronal','sagittal'
ip.addParameter('fontSize',10); % font size for labels

%%% parser
ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end

%%% projections
if strcmp(projection,'sagittal')
    viewAngle = [- 90 0];
elseif strcmp(projection,'coronal')
    viewAngle = [0 0];
else % axial default
    viewAngle = [0 90];
end

%%% fib
DV = fib.fiber_bottom_DV;
AP = fib.fiber_bottom_AP;
ML = fib.fiber_bottom_ML;
roiNums = fib.ROI;

%%% now make figure
if newFig == 1
    figure
end
% colormap
if isnumeric(colormapOption) % if a colormap matrix has been supplied
    cmapping = linspace(min(cmapBounds),max(cmapBounds),size(colormapOption,1));
    cmapcolors = colormapOption;
else
    cmapping = linspace(min(cmapBounds),max(cmapBounds),colormapBins);
    eval(['cmapcolors = ' colormapOption '(' num2str(colormapBins) ');']);
    if strmatch(colormapOption,'redblue') % symmetric colormap for redblue
        cmapping = linspace(-max(abs(cmapBounds)), max(abs(cmapBounds)),colormapBins);
    end
end
if numel(outlineWidth)<size(outlineColor,1)
    outlineWidth = repmat(outlineWidth,1,size(outlineColor,1));
end
% roi order
roiNums = roiNums(:);
roiOrder = [blankBubble(:); setdiff(roiNums,blankBubble)];
for r = 1:numel(roiOrder)
    roinum = roiOrder(r);
    x = ML(roinum);
    y = AP(roinum);
    lineWidth = outlineWidth(outlineWidthID(roinum));
    thisOutlineAlpha = outlineAlpha(outlineAlphaID(roinum));
    thisColorAlpha = colorAlpha(colorAlphaID(roinum));
    
    if ismember(roinum,blankBubble)                            
        thisColor = blankBubbleColor;
    else
        thisVal = values(roinum);
        if thisVal<min(cmapping)
            thisColor = cmapcolors(1,:);
        else
            thisColor = cmapcolors(find(cmapping<=thisVal,1,'last'),:);                
        end
    end
    if lineWidth == -1 
        lineColor = thisColor;
        lineWidth = 1;
    else       
        lineColor = outlineColor(outlineID(roinum),:);  
    end
    depth = DV(roinum);
    scatter3(x,y,depth,bubbleSize,'ok','MarkerEdgeColor',lineColor,'LineWidth',lineWidth,'MarkerFaceColor', thisColor,'MarkerFaceAlpha',thisColorAlpha,'MarkerEdgeAlpha',thisOutlineAlpha);       
    hold on
end
%%% str outline
if ~isempty(str)
    str_outline = get_mask_projection_outlines(str.striatum_mask,str,...
        'apply_str_mask',0,'proj_orientations',{projection});
    if strcmp(projection,'coronal')
        x = str_outline.(projection){1};
        y = ones(size(str_outline.(projection){1}))*mean(AP);
        z = str_outline.(projection){2};
    elseif strcmp(projection,'sagittal')
        x = ones(size(str_outline.(projection){1}))*mean(ML);
        y = str_outline.(projection){1};
        z = str_outline.(projection){2};
    else % axial default
        x = str_outline.(projection){1};
        y = str_outline.(projection){2};
        z = ones(size(str_outline.(projection){1}))*mean(DV);
    end
    plot3(x,y,z,'-k') 
end

set(gca,'FontSize',fontSize)
view(gca,viewAngle)
xlabel('ML')
ylabel('AP')
zlabel('DV')
if includeColorbar == 1
    cb = colorbar(gca);
    colormap(gca,cmapcolors)    
    ticks = get(cb,'ticks');
    ticks = linspace(ticks(1),ticks(end),colorbarTickDec);
    tickLabels = round(linspace(cmapping(1),cmapping(end),colorbarTickDec),3);
    set(cb,'Ticks',ticks,'TickLabels',tickLabels,'FontSize',fontSize,'Location',colorbarLoc)
    if ~isempty(colorbarLabel)
        ylabel(cb,colorbarLabel,'Rotation',270,'FontSize',fontSize)
    end
end




% 
grid off 
axis equal
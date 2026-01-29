% plot_smooth_maps
%
% plot axial and sagittal mean projections with striatal outlines
% str = the returned output from get_striatum_vol_mask.m 
% other_vals = a cell array of binary volumes, from which to draw contours
%
% will return 2 plots
% one symmetric redblue one
% one parula one
function plot_smooth_maps(this_map,str,varargin)

    %%%  parse optional inputs %%%
    ip = inputParser;
    % directories and atlas: see https://github.com/HoweLab/MultifiberLocalization    
    ip.addParameter('outlines',[]); % outlines to plot (see get_mask_projection_outlines)
    ip.addParameter('cmap_bounds',[]);   
    ip.addParameter('cmap_option','redblue');   
    ip.parse(varargin{:});
    for j=fields(ip.Results)'
        eval([j{1} '=ip.Results.' j{1} ';']);
    end
    
    % calculate projections
    axial_vol = get_volume_projection(this_map,'axial','mask',str.striatum_mask);
    sagittal_vol = get_volume_projection(this_map,'sagittal','mask',str.striatum_mask);

    % cmap_bounds: default is symmetric, thresholded at 99.5th prctile of
    % abs(axial_map), redblue map
    if isempty(cmap_bounds)
        cmap_bounds = [-1 1]*prctile(abs(axial_vol(:)),99.5);
    end
        
    % get outlines of axial and sagittal mask projections and add to the 
    % top of the list of outlines to plot
    if isempty(outlines)
        outlines.axial = [];
        outlines.sagittal = [];
    end
    proj_orientations = {'axial','sagittal'};
    str_outlines = get_mask_projection_outlines(str.striatum_mask,str,...
        'apply_str_mask',0,'proj_orientations',proj_orientations);
    for p = 1:numel(proj_orientations)
        outlines.(proj_orientations{p}) = ...
            [str_outlines.(proj_orientations{p});...
            outlines.(proj_orientations{p})];
    end
         
    %%% now plot
    figure('Position',[100 100 450 800])
    
    %%% axial
    subplot(2,1,1)
    hi = imagesc(axial_vol,'XData',str.info.ML,'YData',str.info.AP);
    set(hi, 'AlphaData', ~isnan(axial_vol))
    hold on
    % plot contours (will be black - can change later in illustrator)
    for v = 1:size(outlines.axial,1)
        plot(outlines.axial{v,1},outlines.axial{v,2},'-k') 
    end
    caxis(cmap_bounds)
    colormap(cmap_option)
    set(gca,'XDir','Normal','YDir','normal')
    xlabel('Medial \leftrightarrow Lateral')
    ylabel('Posterior \leftrightarrow Anterior')
    axis equal
    colorbar
    title('Axial')


    %%% sagittal
    subplot(2,1,2)
    hi = imagesc(sagittal_vol,'XData',str.info.AP,'YData',str.info.DV);
    set(hi, 'AlphaData', ~isnan(sagittal_vol))
    hold on
    % plot contours (will be black - can change later in illustrator)
    for v = 1:size(outlines.sagittal,1)
        plot(outlines.sagittal{v,1},outlines.sagittal{v,2},'-k') 
    end
    caxis(cmap_bounds)
    colormap(cmap_option)
    set(gca,'XDir','reverse','YDir','normal') % for whatever reason, I like anterior on the left
    xlabel('Anterior \leftrightarrow Posterior')
    ylabel('Ventral \leftrightarrow Dorsal')
    axis equal  
end    
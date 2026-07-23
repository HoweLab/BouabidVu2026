% plot_contour_maps
%
% Mai-Anh Vu, 2026

function plot_contour_maps(this_map,str,contour_levels,peak_sign,varargin)

    %%%  parse optional inputs %%%
    ip = inputParser;
    % directories and atlas: see https://github.com/HoweLab/MultifiberLocalization    
    ip.addParameter('cmap_option','gray');
    ip.addParameter('peak_color',[1 0 0]);  % leave empty to omit peak
    ip.addParameter('str_color',[1 1 .9]);   
    ip.parse(varargin{:});
    for j=fields(ip.Results)'
        eval([j{1} '=ip.Results.' j{1} ';']);
    end
    
    % get peak
    if peak_sign == 1
        [~,peak_ind] = nanmax(this_map(:));
    elseif peak_sign == -1
        [~,peak_ind] = nanmin(this_map(:));
    end
    [i,j,k] = ind2sub(size(this_map),peak_ind);    
    peak_val = this_map(peak_ind);
    peak_dv = str.info.DV(i);
    peak_ml = str.info.ML(j);
    peak_ap = str.info.AP(k);
    
    % get contour colors
    eval(['contour_colors = ' cmap_option '(ceil(numel(contour_levels)/2));'])
    if rem(numel(contour_levels),2)==0
        contour_colors = [contour_colors; flipud(contour_colors)];
    else
        contour_colors = [contour_colors; flipud(contour_colors(1:end-1,:))];
    end

    % 3D contour volumes
    c_contour = cell(numel(contour_levels),1);
    for i = 1:numel(contour_levels)
        if peak_sign == 1
            c_contour{i} = this_map >= contour_levels(i);
        elseif peak_sign == -1
            c_contour{i} = this_map <= contour_levels(i);
        end
    end
    
    % flip plot order if peak is negative
    if peak_sign == -1
        c_contour = flipud(c_contour);
        contour_colors = flipud(contour_colors);
    end

    % striatal outlines
    proj_orientations = {'axial','sagittal'};
    str_outlines = get_mask_projection_outlines(str.striatum_mask,str,...
        'apply_str_mask',0,'proj_orientations',proj_orientations);
         
    %%% now plot
    figure('Position',[100 100 450 800])
    
    %%% axial
    subplot(2,1,1)
    hold on
    % fill in striatum
    fill(str_outlines.axial{1},str_outlines.axial{2},str_color)
    % fill in contours
    for c = 1:numel(c_contour)
        this_contour = c_contour{c};
        if sum(this_contour(:))>0
            contour_proj = get_volume_projection(this_contour,'axial',...
                'proj','max','mask',str.striatum_mask,'mask_replace',0);                    
            contour_outline = bwboundaries(contour_proj,'noholes');
            for b = 1:numel(contour_outline)
                xy = contour_outline{b};
                x = str.info.ML(xy(:,2));
                y = str.info.AP(xy(:,1));          
                plot(x,y)
                if numel(x)>2
                    fill(x,y,contour_colors(c,:),'FaceAlpha',1,'EdgeColor',[0 0 0])
                end
            end
        end
    end
    % plot peak
    if ~isempty(peak_color)
        plot(peak_ml,peak_ap,'.','MarkerSize',15,'Color',peak_color)
    end
    set(gca,'XDir','Normal','YDir','normal')
    xlabel('Medial \leftrightarrow Lateral')
    ylabel('Posterior \leftrightarrow Anterior')
    colormap(gca,contour_colors)            
    cb = colorbar;
    set(cb,'Ticks',[0 .5 1],'TickLabels',[min(contour_levels) 0 max(contour_levels)])



    %%% sagittal
    subplot(2,1,2)
    hold on
    % fill in striatum
    fill(str_outlines.sagittal{1},str_outlines.sagittal{2},str_color)
    % fill in contours
    for c = 1:numel(c_contour)
        this_contour = c_contour{c};
        if sum(this_contour(:))>0
            contour_proj = get_volume_projection(this_contour,'sagittal',...
                'proj','max','mask',str.striatum_mask,'mask_replace',0);                    
            contour_outline = bwboundaries(contour_proj,'noholes');
            for b = 1:numel(contour_outline)
                xy = contour_outline{b};
                x = str.info.AP(xy(:,2));
                y = str.info.DV(xy(:,1));          
                plot(x,y)
                if numel(x)>2
                    fill(x,y,contour_colors(c,:),'FaceAlpha',1,'EdgeColor',[0 0 0])
                end
            end
        end
    end
    % plot peak
    if ~isempty(peak_color)
        plot(peak_ap,peak_dv,'.','MarkerSize',15,'Color',peak_color)
    end
    set(gca,'XDir','reverse','YDir','normal') % for whatever reason, I like anterior on the left
    xlabel('Anterior \leftrightarrow Posterior')
    ylabel('Ventral \leftrightarrow Dorsal')
    axis equal  
end    
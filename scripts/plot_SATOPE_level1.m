function plot_SATOPE_level1(h, eventdata,fn,fp,Ainc,varargin)

% Function to be called as CallBack from plot_SATOPE_level4
% in order to plot the same as level_4 but for the high-res
% brightness temperatures CMEM ouput level 1.
% This assumes that the level_1 output files are located at
% the same directory as the level_4 and with the classical
% file name convention as 'CMEM'
% 
% * TO DO: implement this function as optional stand-alone.
% (c) 2017 P. Saavedra Garfias, UNIVERSITY OF BONN, GERMANY
% Email: pablosaa@uni-bonn.de
% See LICENSE.TXT
% --------------------------------------------------------
    
    [ho, fo] = gcbo;
    if isnumeric(fo), tt=fo; else tt=fo.Number; end;
    lstrtime = 'UTC Time: %4.1f hr';
    optvars = {'LONLIM','LATLIM','CLIM_h','CLIM_v'};
    for i=1:length(optvars),
        eval([optvars{i} '=[];']);
    end
    for i=1:nargin-5,
        eval([optvars{i} '=varargin{i};']);
    end
    last_usc = max(strfind(fn{tt},'_'));
    Ninc = length(Ainc{tt});
    incang = Ainc{tt};
    for i=1:Ninc,
        fname = sprintf('%s%s%s%02d.nc',fp,'out_level1',...
                        fn{tt}(11:last_usc),incang(i));
        if exist(fname,'file'),
            disp(['Reading ' fname]);
            if i==1,                
                time = ncread(fname,'TIME')/60/60;
                lon = ncread(fname,'LONGITUDE');
                lat = ncread(fname,'LATITUDE');
                if isempty(LONLIM) | isempty(LATLIM),
                    LONLIM = [0.999*min(lon) 1.005*max(lon)];
                    LATLIM = [0.999*min(lat) 1.005*max(lat)];
                end
            end
            tb_h(:,:,:,i) = ncread(fname,'TBH');
            tb_v(:,:,:,i) = ncread(fname,'TBV');            
        else
            warning([':O ' fname ', file not found!']);
            tb_h(:,:,:,i) = NaN;
            tb_v(:,:,:,i) = NaN;
        end
    end
    Ntim = length(time);
    figure(10+tt);
    clf;
    set(gcf,'Position',[169 21 355*Ninc 930]);
    k=1;
    for j=1:Ninc,
        ax1(j,1) = subplot(2,Ninc,j); % 2*j-1
        pixs(j,1) = pcolor(lon,lat,tb_h(:,:,k,j)');
        shading flat;
        title(sprintf('H-pol: \\theta_{inc} = %3.1f deg',incang(j)),...
              'FontSize',15);
        ax1(j,2) = subplot(2,Ninc,j+Ninc); % 2*j
        pixs(j,2) = pcolor(lon,lat,tb_v(:,:,k,j)');
        shading flat;
        title(sprintf('V-pol: \\theta_{inc} = %3.1f deg',incang(j)),...
              'FontSize',15);
        % checking limits for the color-scale:
        if isempty(CLIM_h),
            CLIM_h = [min(min(min(tb_h(:,:,k,j)))),...
                      max(max(max(tb_h(:,:,k,j))))];
        end
        if isempty(CLIM_v),
            CLIM_v = [min(min(min(tb_v(:,:,k,j)))),...
                      max(max(max(tb_v(:,:,k,j))))];
        end
        set(ax1(j,1),'CLim',CLIM_h);
        set(ax1(j,2),'CLim',CLIM_v);
    end  % end loop over inc angles
        
    % Setting subplots properties and colorbar:
    set(ax1,'FontSize',15,'TickDir','out','Box','on',...
            'XLim',LONLIM,'YLim',LATLIM,...
            'XGrid','on','YGrid','on');
    tmp = get(ax1(j,end),'Position');
    hbar = colorbar('SouthOutside');
    set(ax1(j,end),'Position',tmp);  % returning to original
                                    % position
    set(hbar,'Position',[0.28 0.065 0.49 0.012],'FontSize',15);
    set(get(hbar,'Xlabel'),'FontSize',16,'String',...
                      ['L-band Brightness Temperature [K]']);

    % Inserting the time slide GUI
    % Text and Slider for Time variable:
    tti = uicontrol(gcf,'Style','text','Position',[10 910 110 20],...
                    'HorizontalAlignment','Left',...
                    'String',[sprintf(lstrtime,time(k))]);
    if Ntim>1,
        hti = uicontrol(gcf,'Style','slider','Min',1,...
                        'Max',Ntim,'Position',[120 910 150 20],...
                        'SliderStep',[1 1],...
                        'Value',k,'CallBack',...
                        {@ChangeTimeHr,tti,Ninc,pixs,tb_h,tb_v,time});
    end
    return;
end  % function plot_SATOPE_level1
    
function ChangeTimeHr(h,eventdata,tti,Ninc,pixs,tb_h,tb_v,time)
    lstrtime = 'UTC Time: %4.1f hr';
    [ho,fo] = gcbo;
    k=round(get(ho,'Value'));
      if isnumeric(fo), tt=fo; else tt=fo.Number; end;
      for p=1:Ninc,
          set(pixs(p,1),'CData',tb_h(:,:,k,p)');
          set(pixs(p,2),'CData',tb_v(:,:,k,p)');
          set(tti,'String',sprintf(lstrtime,time(k)));
      end
end % function ChangeTimeHr

    
    % end of script

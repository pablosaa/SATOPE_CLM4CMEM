function [ax, pixs] = plot_SATOPE_level4(varargin)
% Script to plot the output of the CLM4CMEM level 4 NetCDF data.
% Multiple files selection possible, every selected file will be 
% plot in a separate Figure Window.
% When more than 3 Incidence angles are present, then only 3 are
% showed in every figure.
% 
% OCTAVE/MATLAB Friendly Script. MATLAB only new versions (R2014b>)
%
% (c) 2017 P. Saavedra, UNIVERSITY OF BONN, Germany 
% Email: pablosaa@uni-bonn.de
% See LICENCE.TXT
% ---------------------------------------------------------------

global ax pixs;
%if nargin==0,
[fname, fpath] = uigetfile({'*.nc','NetCDF'},...  % filter *.nc
                    'CMEM Level_4 NetCDF output data to load?',...
                    '../../SATOPE_CLM4CMEM/output/',...  % default dir
                    'MultiSelect','on');   % Multiple selection allowed
%end
if iscell(fname),
    Nfile = length(fname);
elseif ischar(fname),
    Nfile = 1;
    fname = {sprintf('%s',fname)};
else
    error('No file selected. End!');
    return;
end
if exist('OCTAVE_VERSION','builtin')>0,
    disp('Loading NetCDF packages for Octave...');
    pkg load netcdf
	pixsi = 13;
	if all(OCTAVE_VERSION=='4.2.2'),
		pixsi = pixsi^2;
	end
else
    pixsi = 70;
end


lstrtime = 'UTC Time: %4.1f hr';
for i=1:Nfile,
    disp(['File to open: ', fname{i}]);
    % Loading data:
    le4info{i} = ncinfo([fpath fname{i}]);
    temp = squeeze(struct2cell(le4info{i}.Attributes));
    infomsg{i} = strcat(upper(temp(1,:)),':= ', temp(2,:));
    %npix = ncread([fpath fname{i}],'NPIXEL');
    incang{i} = ncread([fpath fname{i}],'INC_ANGLE');
    Ninc{i} = length(incang{i});
    time{i} = ncread([fpath fname{i}],'TIME')/60/60;  % time in hr UTC
    Ntim{i} = length(time{i});
    %disp(sprintf('Plotting number of pixels: %d\n', npix));
    lon{i} = ncread([fpath fname{i}],'LONGITUDE');
    lat{i} = ncread([fpath fname{i}],'LATITUDE');
    tb_h{i} = ncread([fpath fname{i}],'TBSAT_H');
    tb_h{i}(tb_h{i}<=0) = NaN;
    tb_v{i} = ncread([fpath fname{i}],'TBSAT_V');
    tb_v{i}(tb_v{i}<=0) = NaN;
    % checking limits for the color-scale:
    CLIM_h = [min(min(min(tb_h{i}))), max(max(max(tb_h{i})))];
    CLIM_v = [min(min(min(tb_v{i}))), max(max(max(tb_v{i})))];
    LONLIM = [0.99*min(lon{i}) 1.01*max(lon{i})];
    LATLIM = [0.999*min(lat{i}) 1.005*max(lat{i})];
    % Plotting data:
    k=1;  % for the first time index
    figure(i);
    clf;
    set(gcf,'Position',[169 21 355*Ninc{i} 930]);
    for j=1:Ninc{i},
        ax{i}(j,1) = subplot(2,Ninc{i},j); % 2*j-1
        pixs{i}(j,1) = scatter(lon{i},lat{i},pixsi,tb_h{i}(:,j,k),'o','filled');
        title(sprintf('H-pol: \\theta_{inc} = %3.1f deg',incang{i}(j)),...
              'FontSize',15);
        ax{i}(j,2) = subplot(2,Ninc{i},j+Ninc{i}); % 2*j
        pixs{i}(j,2) = scatter(lon{i},lat{i},pixsi,tb_v{i}(:,j,k),'o','filled');
        title(sprintf('V-pol: \\theta_{inc} = %3.1f deg',incang{i}(j)),...
              'FontSize',15);
    end

    % Text for the Satellite Name:
    tih = uicontrol(gcf,'Style','Text','Position',[400 910 300 30],...
                    'HorizontalAlignment','center',...
                    'String',[le4info{i}.Attributes(3).Name ': ' le4info{i}.Attributes(3).Value],...
                    'FontSize',20);
    % Text for the Information of the souce file:
    tmsg = uicontrol(gcf,'Style','pushbutton','Position',[800 910 80 30],...
                     'String','INFO','CallBack',...
                     {@showtheinfo,infomsg});

    % Button for generate equivalent plot with level_1 data:
    tbu = uicontrol(gcf,'Style','pushbutton','String','Plot Level 1',...
                    'Position',[890 910 100 30],...
                    'CallBack',{@plot_SATOPE_level1,fname,fpath,incang,...
                        LONLIM,LATLIM,CLIM_h,CLIM_v});
    % Text and Slider for Time variable:
    tti{i} = uicontrol(gcf,'Style','text','Position',[10 910 110 20],...
                       'HorizontalAlignment','Left',...
                       'String',[sprintf(lstrtime,time{i}(k))]);
    if Ntim{i}>1,
        hti{i} = uicontrol(gcf,'Style','slider','Min',1,...
                           'Max',Ntim{i},'Position',[120 910 150 20],...
                           'SliderStep',[1 1],...
                           'Value',k,'CallBack',...
                           {@ChangeTimeHr,tti,Ninc,pixs,tb_h,tb_v,time});
    end
    % Text and Input for Pixel Size (in case too small/big in plot):
    hpt = uicontrol(gcf,'Style','Text','String','Pixel Size:',...
                    'Position',[325 910 70 20]);
    hpz = uicontrol(gcf,'Style','Edit','String',num2str(pixsi),...
                    'Position',[400 910 40 20],'CallBack',...
                    {@ChangePixelSize,pixs});

    % Setting subplots properties and colorbar:
    set(ax{i},'FontSize',15,'TickDir','out','Box','on',...
              'XLim',LONLIM,'YLim',LATLIM,...
              'XGrid','on','YGrid','on');
    set(ax{i}(:,1),'CLim',CLIM_h);
    set(ax{i}(:,2),'CLim',CLIM_v);
    tmp = get(ax{i}(j,end),'Position');
    hbar = colorbar('SouthOutside');
    set(ax{i}(j,end),'Position',tmp);  % returning to original
                                    % position
    set(hbar,'Position',[0.28 0.065 0.49 0.012],'FontSize',15);
    set(get(hbar,'Xlabel'),'FontSize',16,'String',...
                      ['L-band Brightness Temperature [K]']);
end   % end over files

    % end function main()
end

function ChangeTimeHr(h,eventdata,tti,Ninc,pixs,tb_h,tb_v,time)
    lstrtime = 'UTC Time: %4.1f hr';
    [ho,fo] = gcbo;
    k=round(get(ho,'Value'));
    if isnumeric(fo), tt=fo; else tt=fo.Number; end;
    for p=1:Ninc{tt},
        set(pixs{tt}(p,1),'CData',tb_h{tt}(:,p,k));
        set(pixs{tt}(p,2),'CData',tb_v{tt}(:,p,k));
        set(tti{tt},'String',sprintf(lstrtime,time{tt}(k)));
    end
    % end function ChangeTimeHr()
end
        
function ChangePixelSize(h,eventdata,pixs)
    [ho,fo] = gcbo;
    pixsi=str2num(get(ho,'String'));
    if isnumeric(fo), tt=fo; else tt=fo.Number; end;
    set(pixs{tt},'SizeData',pixsi);
    % end function ChangePixelSize()
end

function showtheinfo(h,eventdata,infomsg)
    [ho,fo] = gcbo;
    if isnumeric(fo), tt=fo; else tt=fo.Number; end;
    kk= msgbox(infomsg{tt},'Data Information','help');        
    % end function showtheinfo()
end
% end of script

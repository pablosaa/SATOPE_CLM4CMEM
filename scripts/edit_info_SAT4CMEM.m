function edit_info_SAT4CMEM;
% Script to run a GUI which helps to edit/create NetCDF files
% containing information about specific SATELLITE sensors to be
% used as parameters for the satellite operator and CMEM.
%
% (c) 2016 P. Saavedra Garfias, UNIVERSITY OF BONN
% Email: pablosaa@uni-bonn.de
% See: LICENSE.TXT
% ---------------------------------------------------------------

figure(1);
set(1,'Menu','none','Position',[169 21 500 700]);
% selecting NetCDF SATELLITE info file as default?

myuser = getenv('USER');
myhost = getenv('HOSTNAME');
% default SAT
SAT.NPIX = NaN;
SAT.name{1} = 'SATELLITE_name';
SAT.value{1} = 'GPM';
SAT.attribute{1} = {'global','string',''};
SAT.name{2} = 'Orbit_altitude_km';
SAT.value{2} = 407;
SAT.attribute{2} = {'global','number','km'};
SAT.name{3} = 'Orbit_azimuth_deg';
SAT.value{3} = 172;
SAT.attribute{3} = {'global','number','deg'};
SAT.name{4} = 'SENSOR_antenna_m';
SAT.value{4} = 1.5;
SAT.attribute{4} = {'global','number','m'};
SAT.name{5} = 'SENSOR_wavelength_m';
SAT.value{5} = 0.03;
SAT.attribute{5} = {'global','number','m'};
SAT.name{6} = 'THETA_INC';
SAT.value{6} = [30, 40, 50];
SAT.attribute{6} = {'var','number','deg'};
SAT.name{7} = 'LONG';
SAT.value{7} = [];
SAT.attribute{7} = {'var','number','deg'};
SAT.name{8} = 'LATI';
SAT.value{8} = [];
SAT.attribute{8} = {'var','number','deg'};
SAT.name{9} = 'INCLI';
SAT.value{9} = [];
SAT.attribute{9} = {'var','number','deg'};
SAT.name{10} = 'creator contact';
SAT.value{10} = [myuser '@' myhost];
SAT.attribute{10} = {'global','string',''};

axim = axes('Position',[.05 .8 .35 .19]);
logo = imread('TB_SMAPcube.png');
imshow(logo);

about = uicontrol(gcf,'Style','text','Position',[230 530 250 120],...
                  'String',{upper('Satellite Operator Parameters'),...
                    'This is used to create a satellite operator',...
                   'The structure variable SAT contains:'},...
                  'FontSize',10,'FontName','utopia',...
                  'ForegroundColor','b','BackgroundColor','g');

attfr = uipanel('Title','Global Attributes','FontSize',12,...
                'BackgroundColor','w','Position',[0.05 0.45 .9 .3]);
% Name field:
globalatt(1,1) = uicontrol('Parent',attfr,'Style','text','String','Name: ',...
                           'Position',[5 155 180 20],...
                           'HorizontalAlignment','left');

gloablatt(1,2) = uicontrol('Parent',attfr,'Style','edit','String',...
                           'GPM X-band','Tag','SATELLITE_name',...
                           'Position',[190 155 190 20]);
% Orbit altitute field:
globalatt(2,1) = uicontrol('Parent',attfr,'Style','text',...
                           'String','Orbit H [km]: ',...
                           'Position',[5 125 180 20],...
                           'HorizontalAlignment','left');

gloablatt(2,2) = uicontrol('Parent',attfr,'Style','edit','String','407.0',...
                           'Tag','Orbit_altitude_km',...
                           'Position',[190 125 190 20]);
% Orbit azimuth:
globalatt(3,1) = uicontrol('Parent',attfr,'Style','text','String','Orbit Azimuth [deg]: ',...
                           'Position',[5 95 180 20],...
                           'HorizontalAlignment','left');

gloablatt(3,2) = uicontrol('Parent',attfr,'Style','edit','String','172',...
                           'Tag','Orbit_azimuth_deg',...
                           'Position',[190 95 190 20]);
% Sensor antenna:
globalatt(4,1) = uicontrol('Parent',attfr,'Style','text','String','Antenna diameter [m]: ',...
                           'Position',[5 65 180 20],...
                           'HorizontalAlignment','left');

gloablatt(4,2) = uicontrol('Parent',attfr,'Style','edit','String','1.5',...
                           'Tag','SENSOR_antenna_m',...
                           'Position',[190 65 190 20]);
% Sensor wavelength:
globalatt(5,1) = uicontrol('Parent',attfr,'Style','text','String','Sensor wavelength [m]: ',...
                           'Position',[5 35 180 20],...
                           'HorizontalAlignment','left');

gloablatt(5,2) = uicontrol('Parent',attfr,'Style','edit','String','0.03',...
                           'Tag','SENSOR_wavelength_m',...
                           'Position',[190 35 190 20]);
% Optional attribute: contact name
globalatt(6,1) = uicontrol('Parent',attfr,'Style','text','String','Creator contact: ',...
                           'Position',[5 5 180 20],...
                           'HorizontalAlignment','left');

gloablatt(6,2) = uicontrol('Parent',attfr,'Style','edit','String',[myuser '@' myhost],...
                           'Tag','creator contact',...
                           'Position',[190 5 190 20]);
% ----------------
varfr = uipanel('Title','Variables','FontSize',12,...
                'BackgroundColor','w','Position',[0.05 0.18 .9 .25]);

gloablatt(7,1) = uicontrol('Parent',varfr,'Style','text','String',...
                           'Incidence Angles [deg]', 'Position',...
                           [5 90 180 20],...
                           'HorizontalAlignment','left');
gloablatt(7,2) = uicontrol('Parent',varfr,'Style','edit',...
                           'String','[30, 40, 50]',...
                           'Tag','THETA_INC',...
                           'Position',[190 90 200 20]);
gloablatt(8,1) = uicontrol('Parent',varfr,'Style','text','String',...
                           'Longitude +E [deg]', 'Position',...
                           [5 65 180 20],'HorizontalAlignment','left');
gloablatt(8,2) = uicontrol('Parent',varfr,'Style','edit',...
                           'String',' ',...
                           'Tag','LONG',...
                           'Position',[190 65 200 20]);
gloablatt(9,1) = uicontrol('Parent',varfr,'Style','text','String',...
                           'Latitude +N [deg]', 'Position',...
                           [5 40 180 20],'HorizontalAlignment','left');
gloablatt(9,2) = uicontrol('Parent',varfr,'Style','edit',...
                           'String',' ',...
                           'Tag','LATI',...
                           'Position',[190 40 200 20]);
gloablatt(10,1) = uicontrol('Parent',varfr,'Style','text','String',...
                           'Inclination CCW [deg]', 'Position',...
                           [5 15 180 20],'HorizontalAlignment','left');
gloablatt(10,2) = uicontrol('Parent',varfr,'Style','edit',...
                           'String',' ',...
                           'Tag','INCLI',...
                           'Position',[190 15 200 20]);
% ----------------
% Showing Number of Pixels and Angles:

% Button to load existing NetCDF SAT file:

btload = uicontrol(gcf,'Style','pushbutton','String','Load SAT file',...
                   'Position',[30 50 100 30],'CallBack',...
                   {@load_file,gloablatt},'Tag','Load');

btsave = uicontrol(gcf,'Style','pushbutton','String','Save SAT file',...
                   'Position',[150 50 100 30],'CallBack',{@save_file, ...
                    gloablatt,SAT},'Tag','Save');

btmove = uicontrol(gcf,'Style','pushbutton','String','to WorkSpace',...
                   'Position',[270 50 100 30],'CallBack',{@save_file, ...
                    gloablatt,SAT},'Tag','Move');

btexit = uicontrol(gcf,'Style','pushbutton','String','Exit',...
                   'Position',[390 50 100 30],'CallBack','close(gcf)');
end
% end of script

% Function to up-load old parameters stored in a NetCDf file:
function load_file(aa,bb,gloablatt)
    [ho, fo] = gcbo;
    if isnumeric(fo), tt=fo; else tt=fo.Number; end
    SAT = retrieve_SATinputdata;
    nn = length(gloablatt); %numfields(SAT);
    %namef = fieldnames(SAT);
    for i=1:nn,
        namex = get(gloablatt(i,2),'Tag');
        tagnum = find(strcmp(SAT.name,namex));
        value = SAT.value{tagnum}; %getfield(SAT,namef{i+5});

        switch SAT.attribute{tagnum}{1},
          case 'var',            
            if length(value)<8,
                namex = ['[' sprintf('%5.1f',value) ']'];
            else
                namex = sprintf('SAT.value{%d}',tagnum);
            end           
          case 'global',
            if isnumeric(value),
                namex = num2str(value);
            else
                namex = value;
            end
        end
        set(gloablatt(i,2),'String',namex);
    end
    assignin('base','SAT',SAT);
    % TODO: include NPIX and NTHETA?
end

% Function to store the actual parameters in GUI:
function save_file(aa,bb,gloablatt,SAT)
    
    MYOCT = exist('OCTAVE_VERSION','builtin');
    [ho, fo] = gcbo;
    if isnumeric(fo), tt=fo; else tt=fo.Number; end

    for i=1:length(gloablatt),
        value = strtrim(get(gloablatt(i,2),'String'));
        namex = get(gloablatt(i,2),'Tag');
        idx = find(strcmp(SAT.name,namex));
        ilabel = SAT.name{idx};
        if isempty(value),
            set(gloablatt(i,2),'BackgroundColor','r');
            hh = msgbox(['Please write a valid value for "',...
                         ilabel '" in the field box!'],'TRY AGAIN!',...
                        'error','modal');
            if ~MYOCT,        uiwait(hh);      end
            set(gloablatt(i,2),'BackgroundColor','w');
            return;
        end
        
        switch SAT.attribute{idx}{1},
          case 'var'           
            tmp = evalin('base',value,'NaN;');
            if isnan(tmp),
                set(gloablatt(i,2),'BackgroundColor','r');
                hh = msgbox(['Field ' SAT.name{idx} ' has no assigment or variable ' value ' does not exist!'],'TRY AGAIN!',...
                            'error','modal');
                if ~MYOCT,        uiwait(hh);      end
                set(gloablatt(i,2),'BackgroundColor','w');
                return;
            else
                value = tmp;
            end
          case 'global'
            if strcmp(SAT.attribute{idx}{2},'number'),
                value = str2num(value);
            end
        end
        if isvector(value),
            SAT.value{idx} = value;
        else
            set(gloablatt(i,2),'BackgroundColor','r');
            hh = msgbox(['Field ' SAT.name{idx} ' must be a vector!'],'TRY AGAIN!',...
                        'error','modal');
            if ~MYOCT,        uiwait(hh);      end
            set(gloablatt(i,2),'BackgroundColor','w');
            return;           
        end
       
    end

    idx = find(strcmp(SAT.name,'LONG'));
    Nlon = length(SAT.value{idx});
    idx = find(strcmp(SAT.name,'LATI'));
    Nlat = length(SAT.value{idx});
    idx = find(strcmp(SAT.name,'INCLI'));
    Ninc = length(SAT.value{idx});
    
    if (Nlon ~= Nlat | Nlon ~= Ninc | Nlat ~= Ninc),
        set(gloablatt([8:10],2),'BackgroundColor','y');
        hh = msgbox(['Number of elements for Longitude ('...
                     num2str(Nlon) '), Latitude ('...
                     num2str(Nlat) ') and Inclination ('...
                     num2str(Ninc) ') muss be equal!'],'ERROR!',...
                    'error','modal');
        if ~MYOCT, uiwait(hh); end
        set(gloablatt([8:10],2),'BackgroundColor','w','String','');
    else
        SAT.NPIX = Nlon;
        if strcmp(get(ho,'Tag'),'Save'),
            [fname, fpath] = uiputfile('*.nc','Save NetCDF SAT info','./');
            status = create_SATinput_nc([fpath fname],SAT);
            qq=questdlg({['File ' fname ' saved!'],...
                         'Do you want to put the new structure SAT to workspace too?'},'Copy new SAT to workspace');
        elseif strcmp(get(ho,'Tag'),'Move'),
            % To only move values from GUI to Workspace:
            set(ho,'BackgroundColor','r','FontSize',13);
            qq='Yes';
            disp('Structure variable SAT has been updated in Workspace!');
            set(ho,'BackgroundColor','w','FontSize',10);
        end
        switch qq,
          case 'Yes',
            assignin('base','SAT',SAT);
          case 'No',
          case 'Cancel',
          otherwise,
        end
    end

    return;
end



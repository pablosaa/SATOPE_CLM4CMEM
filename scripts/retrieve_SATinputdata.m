function [SAT,varargout]=retrieve_SATinputdata(varargin)

% Function to retrieve the information from SATELLITE NetCDF file
% used as input for the CLM4CMEM toolbox. This includes all
% necessary information for the Satellite Operator.
% 
% The output structure SAT contains the following fields:
% SAT.NPIX  : Number of Satellite pixels to save
% SAT.name  : String with the name of the sensor e.g. 'SMOS'
% SAT.value : Value corresponding to the SAT.name field
% SAT.attribute: cell array with "name's" attributes like:
%    ={'attribute','type','unit'}, e.g. ={'global','numter','km'}
%
% The structure SAT fields .name, .value and .attribute can have any
% different assignations as cells, followiing the example:
%
% SAT.NPIX = 99;
% SAT.name{1} = 'SATELLITE_name';
% SAT.value{1} = 'GPM';
% SAT.attribute{1} = {'global','string',''};
% SAT.name{2} = 'Orbit_altitude_km';
% SAT.value{2} = 407;
% SAT.attribute{2} = {'global','number','km'};
% SAT.name{3} = 'Orbit_azimuth_deg';
% SAT.value{3} = 172;
% SAT.attribute{3} = {'global','number','deg'};
% SAT.name{4} = 'SENSOR_antenna_m';
% SAT.value{4} = 1.5;
% SAT.attribute{4} = {'global','number','m'};
% SAT.name{5} = 'SENSOR_wavelength_m';
% SAT.value{5} = 0.03;
% SAT.attribute{5} = {'global','number','m'};
% SAT.name{6} = 'THETA_INC';
% SAT.value{6} = [30, 40, 50];
% SAT.attribute{6} = {'var','number','deg'};
% ...
% other names are: 'LONG', 'LATI', 'INCLI', 'creator contact', etc.
% with attributes: 'var' for the first three and 'global' for the last one.
%
% (c) 2016 P. Saavedra Garfias, UNIVERSITY OF BONN
% Email: pablosaa@uni-bonn.de
% See: LICENSE.TXT
% ---------------------------------------------------------------

    if exist('OCTAVE_VERSION','builtin'),
        piv = '_';
        pkg load netcdf;
    else
        piv = '.';
    end
    SAT=[];
    switch nargin,
      case 3,
        error(['wrong number of input arguments!']);
      case 1,
        fname = varargin{1};
        fpath = './';
      case 2,
        fname = varargin{1};
        fpath = varargin{2};
      otherwise,
        [fname,fpath] = uigetfile('*.nc',...
                                  'SATELLITE NetCDF info file',...
                                  '../../SATOPE_CLM4CMEM/forcing');
        if isnumeric(fname),
            warning('File not selected!');
            return;
        end
    end
    if ~exist([fpath fname],'file'),
        error(['File ' fpath fname ' does not exist!']);
    end
    
    disp(fname);
    satfields = {'NPIX','name','value','attribute'};
    
    % Opening the NetCDF file:
    eval(['nc = netcdf' piv 'open([fpath fname],''NOWRITE'');']);
    % getting information from NetCDF file:
    eval(['[NDIMS,NVARS,NGATTS,UNLIM] = netcdf' piv 'inq(nc);']);
    % assigning variables:
    idx = 1;

    for varid=0:NVARS-1,
        eval(['SAT.name{idx} = netcdf' piv 'inqVar(nc,varid);']);
        eval(['valor = netcdf' piv 'getVar(nc,varid);']);
        SAT.value{idx} = double(valor);
        SAT.attribute{idx} = {'var','number','deg'};
        idx = idx + 1;          %idx++;
    end;
    % assigning global attributes:
    eval(['gloid = netcdf' piv 'getConstant(''global'');']);
    for varid=0:NGATTS-1,
        eval(['namex = netcdf' piv 'inqAttName(nc,gloid,varid);']);
        eval(['valor = netcdf' piv 'getAtt(nc,gloid,namex);']);
        SAT.name{idx} = namex;
        if isnumeric(valor),
            SAT.value{idx} = double(valor);
            SAT.attribute{idx} = {'global','number',''};
        else
            SAT.value{idx} = valor;
            SAT.attribute{idx} = {'global','string',''};
        end
        
        idx = idx +1;  %idx++;
    end
    % Reading variables:
    
    eval(['[namex, SAT.NPIX] = netcdf' piv 'inqDim(nc,1);']);
    if(nargout>1),
        eval(['[namex, NTHETA_INC] = netcdf' piv 'inqDim(nc,0);']);
        varargout{1}=NTHETA_INC;
    end
    % closing NetCDF file:
    eval(['netcdf' piv 'close(nc);']);
    
end
% end of script

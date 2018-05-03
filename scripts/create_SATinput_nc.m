function status=create_SATinput_nc(fname,SAT)
% Function to create a SATELLITE NetCDF file as input for the CLM4CMEM 
% toolbox, including all necessary information for the Satellite Operator
% 
% The structure SAT must contain the following fields:
% SAT.NPIX  : Number of Satellite pixels to save
% SAT.name  : String with the name of the sensor e.g. 'SMOS'
% SAT.value : Value corresponding to the SAT.name field
% SAT.attribute: cell array with "name's" attributes like:
%    ={'attribute','type','unit'}, e.g. ={'global','numter','km'}
%
% The structure SAT fields .name, .value and .attribute can have any
% different assignations as cells, followiing the example:
%
% SAT.NPIX = Number of Satellite pixels to save
% SAT.name{1} = 'SATELLITE_name';
% SAT.value{1} = 'GPM';  % e.g. SMOS, SMAP
% SAT.attribute{1} = {'global','string',''};
% SAT.name{2} = 'Orbit_altitude_km';
% SAT.value{2} = 407;  % [km]
% SAT.attribute{2} = {'global','number','km'};
% SAT.name{3} = 'Orbit_azimuth_deg';
% SAT.value{3} = 172;   % [deg] from North clockwise
% SAT.attribute{3} = {'global','number','deg'};
% SAT.name{4} = 'SENSOR_antenna_m';
% SAT.value{4} = 1.5;
% SAT.attribute{4} = {'global','number','m'};
% SAT.name{5} = 'SENSOR_wavelength_m';
% SAT.value{5} = 0.03;
% SAT.attribute{5} = {'global','number','m'};
% SAT.name{6} = 'THETA_INC';
% SAT.value{6} = [30, 40, 50];   % Incidence angle from nadir
% SAT.attribute{6} = {'var','number','deg'};
% ...
% other names are: 'LONG', 'LATI', 'INCLI', 'creator contact', etc.
% attributes are 'var' for the first three and 'global' for the last one.
%
% 'INCLI' is Inclination angle of the pixel's footprint CounterClockwise in deg
%
% (c) 2016 P. Saavedra Garfias, UNIVERSITY OF BONN
% Email: pablosaa@uni-bonn.de
% See: LICENSE.TXT
% ---------------------------------------------------------------

status =0;
if exist('OCTAVE_VERSION','builtin'),
    piv = '_';
    pkg load netcdf;
else
    piv = '.';
end

varnam = {'THETA_INC','LONG','LATI','INCLI'};
varids = {'th_varid','lon_varid','lat_varid','incl_varid'};
gloatt = {'SATELLITE_name','Orbit_altitude_km','Orbit_azimuth_deg',...
         'SENSOR_antenna_m','SENSOR_wavelength_m','creator contact'};

% Creating NetCDF file:
eval(['nc = netcdf' piv 'create(fname,''clobber'');']);

idx = find(strcmp(SAT.name,'THETA_INC'));
NINC = length(SAT.value{idx});
eval(['ninc_dimid = netcdf' piv 'defDim(nc,''NINC'',NINC);'])
eval(['npix_dimid = netcdf' piv 'defDim(nc,''NPIXEL'',SAT.NPIX);'])

eval(['th_varid = netcdf' piv 'defVar(nc,''THETA_INC'',''float'',ninc_dimid);'])
eval(['netcdf' piv 'putAtt(nc,th_varid,''unit'',''degrees_nadir'');'])
eval(['netcdf' piv 'putAtt(nc,th_varid,''long_name'',''Incidence Angle'');'])
eval(['netcdf' piv 'putAtt(nc,th_varid,''short_name'',''theta_inc'');'])

% Longitude:
eval(['lon_varid = netcdf' piv 'defVar(nc,''LONG'',''float'',npix_dimid);'])
eval(['netcdf' piv 'putAtt(nc,lon_varid,''unit'',''degree_east'');'])
eval(['netcdf' piv 'putAtt(nc,lon_varid,''long_name'',''FootPrint Longitude'');'])
eval(['netcdf' piv 'putAtt(nc,lon_varid,''short_name'',''foprt_lon'');'])

% Latitude:
eval(['lat_varid = netcdf' piv 'defVar(nc,''LATI'',''float'',npix_dimid);'])
eval(['netcdf' piv 'putAtt(nc,lat_varid,''unit'',''degree_north'');'])
eval(['netcdf' piv 'putAtt(nc,lat_varid,''long_name'',''FootPrint Latitude'');'])
eval(['netcdf' piv 'putAtt(nc,lat_varid,''short_name'',''foprt_lat'');'])

% inclination
eval(['incl_varid = netcdf' piv 'defVar(nc,''INCLI'',''float'',npix_dimid);'])
eval(['netcdf' piv 'putAtt(nc,incl_varid,''unit'',''degree_equator'');'])
eval(['netcdf' piv 'putAtt(nc,incl_varid,''long_name'',''FootPrint Inclination'');'])
eval(['netcdf' piv 'putAtt(nc,incl_varid,''short_name'',''foprt_incl'');'])

eval(['glo_id = netcdf' piv 'getConstant(''global'');'])
for i=1:length(gloatt),
    idx=find(strcmp(SAT.name,gloatt{i}));
    eval(['netcdf' piv 'putAtt(nc,glo_id,gloatt{i},SAT.value{idx});']);
end
% $$$ eval(['netcdf' piv 'putAtt(nc,glo_id,''Orbit_altitude_km'',SAT.orbit);'])
% $$$ eval(['netcdf' piv 'putAtt(nc,glo_id,''Orbit_azimuth_deg'',SAT.azimuth);'])
% $$$ eval(['netcdf' piv 'putAtt(nc,glo_id,''SENSOR_antenna_m'',SAT.antenna);'])
% $$$ eval(['netcdf' piv 'putAtt(nc,glo_id,''SENSOR_wavelength_m'',SAT.wavelength);'])
% $$$ eval(['netcdf' piv 'putAtt(nc,glo_id,''creator_contact'','pablosaa@uni-bonn.de');'])

eval(['netcdf' piv 'endDef(nc);'])

for i=1:length(varnam),
    idx = find(strcmp(SAT.name,varnam{i}));
    eval(['netcdf' piv 'putVar(nc,' varids{i} ',SAT.value{idx});']);
end
% $$$ eval(['netcdf' piv 'putVar(nc,th_varid,SAT.theta);'])
% $$$ eval(['netcdf' piv 'putVar(nc,lon_varid,SAT.lon);'])
% $$$ eval(['netcdf' piv 'putVar(nc,lat_varid,SAT.lat);'])
% $$$ eval(['netcdf' piv 'putVar(nc,incl_varid,SAT.incl);'])

eval(['netcdf' piv 'close(nc);'])

status=1;

return;
% end of function

% Function to create a SATELLITE NetCDF file as input for the CLM4CMEM 
% toolbox, including all necessary information for the Satellite Operator
% 
% The structure SAT must contain the following fields:
% SAT.name : String with the name of the sensor e.g. 'SMOS'
% SAT.orbit: Orbit altitude of the satellite in km, e.g. 750
% SAT.azimuth: Orbit azimuth in deg, e.g. 172 for decending orbit.
% SAT.antenna: Antenna diameter in meters, e.g. 6 for SMAP alike.
% SAT.wavelength: Wavelength of the sensor in cm, e.g. 0.21
% SAT.theta : Vector of incidence angles in deg e.g. [30 40 50]
% SAT.NPIX  : Number of Satellite pixels to save
% SAT.lon : Longitude of the pixel in deg
% SAT.lat : Latitude of the pixel in deg
% SAT.incl : Inclination angle of the pixel's footprint CCW in deg

% (c) 2016 P. Saavedra Garfias (pablosaa@uni-bonn.de) UNI BONN

function status=create_SATinputdata(fname,SAT)

nc = netcdf_create(fname,'clobber');

NINC = length(SAT.theta);
ninc_dimid = netcdf_defDim(nc,'NINC',NINC);
npix_dimid = netcdf_defDim(nc,'NPIXEL',SAT.NPIX);

th_varid = netcdf_defVar(nc,'THETA_INC','float',ninc_dimid);
netcdf_putAtt(nc,th_varid,'unit','degrees_nadir');
netcdf_putAtt(nc,th_varid,'long_name','Incidence Angle');
netcdf_putAtt(nc,th_varid,'short_name','theta_inc');

% Longitude:
lon_varid = netcdf_defVar(nc,'LONG','float',npix_dimid);
netcdf_putAtt(nc,lon_varid,'unit','degree_east');
netcdf_putAtt(nc,lon_varid,'long_name','FootPrint Longitude');
netcdf_putAtt(nc,lon_varid,'short_name','foprt_lon');

% Latitude:
lat_varid = netcdf_defVar(nc,'LATI','float',npix_dimid);
netcdf_putAtt(nc,lat_varid,'unit','degree_north');
netcdf_putAtt(nc,lat_varid,'long_name','FootPrint Latitude');
netcdf_putAtt(nc,lat_varid,'short_name','foprt_lat');

% inclination
incl_varid = netcdf_defVar(nc,'INCLI','float',npix_dimid);
netcdf_putAtt(nc,incl_varid,'unit','degree_equator');
netcdf_putAtt(nc,incl_varid,'long_name','FootPrint Inclination');
netcdf_putAtt(nc,incl_varid,'short_name','foprt_incl');


glo_id = netcdf_getConstant("global");
netcdf_putAtt(nc,glo_id,'SATELLITE_name',SAT.name);
netcdf_putAtt(nc,glo_id,'Orbit_altitude_km',SAT.orbit);
netcdf_putAtt(nc,glo_id,'Orbit_azimuth_deg',SAT.azimuth);
netcdf_putAtt(nc,glo_id,'SENSOR_antenna_m',SAT.antenna);
netcdf_putAtt(nc,glo_id,'SENSOR_wavelength_m',SAT.wavelength);
netcdf_putAtt(nc,glo_id,'creator_contact','pablosaa@uni-bonn.de');

netcdf_endDef(nc);

netcdf_putVar(nc,th_varid,SAT.theta);
netcdf_putVar(nc,lon_varid,SAT.lon);
netcdf_putVar(nc,lat_varid,SAT.lat);
netcdf_putVar(nc,incl_varid,SAT.incl);

netcdf_close(nc);

status=1;

return;
% end of function

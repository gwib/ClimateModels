%% constants
KtoC = -273.15;
latLim = [58, 84];
lonLim = [-91, 29];
%%
%NCEP reanalysis Temperature obs 1
ncname_NCEP='./data/NCEP_reanalysis_air.2m.gauss.1976_2005_remap_mon.nc';
%Display the metadata for the file (go though to familiarize with the format)
ncdisp(ncname_NCEP);

%Get file information (and show dimension names and variable names)
finfo=ncinfo(ncname_NCEP);
dimNames = {finfo.Dimensions.Name};
varNames = {finfo.Variables.Name};
disp(dimNames);
disp(varNames);

%Load the variables
dlon_NCEP = ncread(ncname_NCEP,'lon');
dlat_NCEP = ncread(ncname_NCEP,'lat');
dtime_NCEP = ncread(ncname_NCEP,'time');
TEMP_obs_NCEP = ncread(ncname_NCEP,'air');
TEMP_obs_NCEP = TEMP_obs_NCEP(:,:,end-311:end); % getting out all temp values since Jan 1980
TEMP_obs_NCEP = TEMP_obs_NCEP + KtoC; % K to C
% in Kelvin, daily avg

for i=1:length(dlon_NCEP)
  dlatx_NCEP(i,:)=dlat_NCEP;
end
for i=1:length(dlat_NCEP)
  dlonx_NCEP(:,i)=dlon_NCEP;
end
%dlatx_NCEP = dlatx_NCEP';
%dlonx_NCEP = dlonx_NCEP';

% Visually check data by plotting the total period mean
NCEP_obs_temp = nanmean(TEMP_obs_NCEP,3); 

figure
hold on
axesm miller
worldmap world
geoshow(dlatx_NCEP,dlonx_NCEP,NCEP_obs_temp,'displaytype','texturemap');
load coast
colorbar
caxis([-30 30])
title('Total period Mean Temperature CCSM4_historical')

figure
hold on
axesm miller
worldmap(latLim, lonLim)
geoshow(dlatx_NCEP,dlonx_NCEP,NCEP_obs_temp,'displaytype','texturemap');
load coastlines
geoshow(coastlat, coastlon)
colorbar
caxis([-30 30])
title('Total period Mean Temperature CCSM4_historical')


%% HadCrut data - obs 2
% ncname_hadCRUT='./data/HadCRUT.4.6.0.0.median_remap.nc';
% %Display the metadata for the file (go though to familiarize with the format)
% ncdisp(ncname_hadCRUT);
% 
% %Get file information (and show dimension names and variable names)
% finfo=ncinfo(ncname_hadCRUT);
% dimNames = {finfo.Dimensions.Name};
% varNames = {finfo.Variables.Name};
% disp(dimNames);
% disp(varNames);
% 
% %Load the variables
% dlon_hadCRUT = ncread(ncname_hadCRUT,'lon');
% dlat_hadCRUT = ncread(ncname_hadCRUT,'lat');
% dtime_hadCRUT = ncread(ncname_hadCRUT,'time');
% TEMP_anomaly_hadCRUT = ncread(ncname_hadCRUT,'temperature_anomaly');

%% Global model data
ncname_tas_Aamon='./data/tas_Amon_EC-EARTH_decadal2000_r3i3p1_198001-201012.nc';
%Display the metadata for the file (go though to familiarize with the format)
ncdisp(ncname_tas_Aamon);

%Get file information (and show dimension names and variable names)
finfo=ncinfo(ncname_tas_Aamon);
dimNames = {finfo.Dimensions.Name};
varNames = {finfo.Variables.Name};
disp(dimNames);
disp(varNames);

%Load the variables
dlon_tas_Aamon = ncread(ncname_tas_Aamon,'lon');
dlat_tas_Aamon = ncread(ncname_tas_Aamon,'lat');
dtime_tas_Aamon = ncread(ncname_tas_Aamon,'time');
dtime = dtime_tas_Aamon(1:312); % global reference time
TEMP_tas_Aamon = ncread(ncname_tas_Aamon,'tas');
TEMP_tas_Aamon = TEMP_tas_Aamon(:,:,1:312);
TEMP_tas_Aamon = TEMP_tas_Aamon + KtoC;
%  monthly avg

for i=1:length(dlon_tas_Aamon)
  dlatx_tas_Aamon(i,:)=dlat_tas_Aamon;
end
for i=1:length(dlat_tas_Aamon)
  dlonx_tas_Aamon(:,i)=dlon_tas_Aamon;
end
%dlatx_NCEP = dlatx_NCEP';
%dlonx_NCEP = dlonx_NCEP';

% Visually check data by plotting the total period mean
TEMP_tas_Aamon_mean = nanmean(TEMP_tas_Aamon,3); 

figure
hold on
axesm miller
worldmap world
geoshow(dlatx_tas_Aamon,dlonx_tas_Aamon,TEMP_tas_Aamon_mean,'displaytype','texturemap');
load coast
colorbar
caxis([-30 30])
title('Total period Mean Temperature Amon')

figure
hold on
axesm miller
worldmap(latLim, lonLim)
geoshow(dlatx_tas_Aamon,dlonx_tas_Aamon,TEMP_tas_Aamon_mean,'displaytype','texturemap');
load coastlines
geoshow(coastlat, coastlon)
colorbar
caxis([-30 30])
title('Total period Mean Temperature Amon')


%% RCM data
ncname_tasArcHistorical='./data/tas_ARC-44_ICHEC-EC-EARTH_historical_r3i1p1_DMI-HIRHAM5_v1_mon_19760101-20051231_remap.nc';
%Display the metadata for the file (go though to familiarize with the format)
ncdisp(ncname_tasArcHistorical);

%Get file information (and show dimension names and variable names)
finfo=ncinfo(ncname_tasArcHistorical);
dimNames = {finfo.Dimensions.Name};
varNames = {finfo.Variables.Name};
disp(dimNames);
disp(varNames);

%Load the variables
dlon_tasArcHistorical = ncread(ncname_tasArcHistorical,'lon');
dlat_tasArcHistorical = ncread(ncname_tasArcHistorical,'lat');
dtime_tasArcHistorical = ncread(ncname_tasArcHistorical,'time');
tas_tasArcHistorical = ncread(ncname_tasArcHistorical,'tas');
tas_tasArcHistorical = tas_tasArcHistorical(:,:,end-311:end);
tas_tasArcHistorical = tas_tasArcHistorical + KtoC;
% Kelvin, monthly(?)

for i=1:length(dlon_tasArcHistorical)
  dlatx_tasArcHistorical(i,:)=dlat_tasArcHistorical;
end
for i=1:length(dlat_tasArcHistorical)
  dlonx_tasArcHistorical(:,i)=dlon_tasArcHistorical;
end
%dlatx_NCEP = dlatx_NCEP';
%dlonx_NCEP = dlonx_NCEP';

% Visually check data by plotting the total period mean
tas_tasArcHistorical_mean = nanmean(tas_tasArcHistorical,3); 

figure
hold on
axesm miller
worldmap world
geoshow(dlatx_tasArcHistorical,dlonx_tasArcHistorical,tas_tasArcHistorical_mean,'displaytype','texturemap');
load coast
colorbar
caxis([-30 30])
title('Total period Mean Temperature Amon')

figure
hold on
axesm miller
worldmap(latLim, lonLim)
geoshow(dlatx_tasArcHistorical,dlonx_tasArcHistorical,tas_tasArcHistorical_mean,'displaytype','texturemap');
load coastlines
geoshow(coastlat, coastlon)
colorbar
caxis([-30 30])
title('Total period Mean Temperature Amon')
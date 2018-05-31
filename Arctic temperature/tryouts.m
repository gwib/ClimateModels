%% CPC data
ncname_CPC='./data/cpc_global_temp.1979_2005_remap.nc';
%Display the metadata for the file (go though to familiarize with the format)
ncdisp(ncname_CPC);

%Get file information (and show dimension names and variable names)
finfo=ncinfo(ncname_CPC);
dimNames = {finfo.Dimensions.Name};
varNames = {finfo.Variables.Name};
disp(dimNames);
disp(varNames);

%Load the variables
dlon_CRU = ncread(ncname_CPC,'lon');
dlat_CRU = ncread(ncname_CPC,'lat');
dtime_CRU = ncread(ncname_CPC,'time');
%TEMP_anomaly_CRU = ncread(ncname_CPC,'temperature_anomaly');


%% HadCrut data
ncname_hadCRUT='./data/HadCRUT.4.6.0.0.median_remap.nc';
%Display the metadata for the file (go though to familiarize with the format)
ncdisp(ncname_hadCRUT);

%Get file information (and show dimension names and variable names)
finfo=ncinfo(ncname_hadCRUT);
dimNames = {finfo.Dimensions.Name};
varNames = {finfo.Variables.Name};
disp(dimNames);
disp(varNames);

%Load the variables
dlon_hadCRUT = ncread(ncname_hadCRUT,'lon');
dlat_hadCRUT = ncread(ncname_hadCRUT,'lat');
dtime_hadCRUT = ncread(ncname_hadCRUT,'time');
TEMP_anomaly_hadCRUT = ncread(ncname_hadCRUT,'temperature_anomaly');

%% HadCrut data
ncname_tasAmon='./data/tas_Amon_EC-EARTH_decadal2000_r3i3p1_198001-201012.nc';
%Display the metadata for the file (go though to familiarize with the format)
ncdisp(ncname_tasAmon);

%Get file information (and show dimension names and variable names)
finfo=ncinfo(ncname_tasAmon);
dimNames = {finfo.Dimensions.Name};
varNames = {finfo.Variables.Name};
disp(dimNames);
disp(varNames);

%Load the variables
dlon_tasAmon = ncread(ncname_tasAmon,'lon');
dlat_tasAmon = ncread(ncname_tasAmon,'lat');
dtime_tasAmon = ncread(ncname_tasAmon,'time');
tas_tasAmon = ncread(ncname_tasAmon,'tas');


%% tas ARC data
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
dlat_tasArcHistorical = ncread(ncname_tasrcHistorical,'lat');
dtime_tasArcHistorical = ncread(ncname_tasArcHistorical,'time');
tas_tasArcHistorical = ncread(ncname_tasArcHistorical,'tas');

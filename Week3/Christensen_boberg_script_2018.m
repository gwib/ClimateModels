clear all

% % TIME INFO

% CCSM4: 185001-200512 - remapped
% EC-earth: 185001-200912 - remapped
% HadGEM2: 185912-200511  - remapped
% CRU (OBS): 196101-200012

ncname='./../data/tas_Amon_CCSM4_historical_r1i1p1_185001-200512_remap_1961_2000.nc';

%Display the metadata for the file (go though to familiarize with the format)
ncdisp(ncname);
%Get file information (and show dimension names and variable names)
finfo=ncinfo(ncname);
dimNames = {finfo.Dimensions.Name};
varNames = {finfo.Variables.Name};
disp(dimNames);
disp(varNames);
%Load the variables
dlon = ncread(ncname,'lon');
dlat = ncread(ncname,'lat');
dtime = ncread(ncname,'time');
CCSM4_historical_temp = ncread(ncname,'tas');

CCSM4_historical_temp = CCSM4_historical_temp - 273.15; % C to K

for i=1:length(dlon)
  dlatx(i,:)=dlat;
end
for i=1:length(dlat)
  dlonx(:,i)=dlon;
end

% Visually check data by plotting the total period mean
CCSM4_historical_temp2 = nanmean(CCSM4_historical_temp,3); 

figure
hold on
axesm miller
worldmap world
geoshow(dlatx,dlonx,CCSM4_historical_temp2,'displaytype','texturemap');
load coast
colorbar
caxis([-30 30])
title('Total period Mean Temperature CCSM4_historical')

%%

ncname='./../data/tas_Amon_EC-EARTH_DMI_historical_r3i1p1_185001-200912_remap_1961_2000.nc'; % insert correct name
%Display the metadata for the file (go though to familiarize with the format)
ncdisp(ncname);
%Get file information (and show dimension names and variable names)
finfo=ncinfo(ncname);
dimNames = {finfo.Dimensions.Name};
varNames = {finfo.Variables.Name};
disp(dimNames);
disp(varNames);
%Load the variables
dlon = ncread(ncname,'lon');
dlat = ncread(ncname,'lat');
dtime = ncread(ncname,'time');
EC_EARTH_DMI_historical_temp = ncread(ncname,'var167');

EC_EARTH_DMI_historical_temp = EC_EARTH_DMI_historical_temp - 273.15; % C to K

% Visually check data by plotting the total period mean
EC_EARTH_DMI_historical_temp2 = nanmean(EC_EARTH_DMI_historical_temp,3); 

figure
hold on
axesm miller
worldmap world
geoshow(dlatx,dlonx,EC_EARTH_DMI_historical_temp2,'displaytype','texturemap');
load coast
colorbar
caxis([-30 30])
title('Total period Mean Temperature EC_EARTH_DMI_historical')



%%

ncname='./../data/tas_Amon_HadGEM2-CC_historical_r1i1p1_185912-200511_remap_1961_2000.nc'; % insert correct name
%Display the metadata for the file (go though to familiarize with the format)
ncdisp(ncname);
%Get file information (and show dimension names and variable names)
finfo=ncinfo(ncname);
dimNames = {finfo.Dimensions.Name};
varNames = {finfo.Variables.Name};
disp(dimNames);
disp(varNames);
%Load the variables
dlon = ncread(ncname,'lon');
dlat = ncread(ncname,'lat');
dtime = ncread(ncname,'time');
HadGEM2_CC_historical_temp = ncread(ncname,'tas');

HadGEM2_CC_historical_temp = HadGEM2_CC_historical_temp - 273.15; % C to K

% Visually check data by plotting the total period mean
HadGEM2_CC_historical_temp2 = nanmean(HadGEM2_CC_historical_temp,3); 

figure
hold on
axesm miller
worldmap world
geoshow(dlatx,dlonx,HadGEM2_CC_historical_temp2,'displaytype','texturemap');
load coast
colorbar
caxis([-30 30])
title('Annual Periodical Mean Temp HadGEM2_CC_historical')

%%

ncname='./../data/cru_ts_3_10.1961.2000.tmp.dat.nc'; % insert correct name
%Display the metadata for the file (go though to familiarize with the format)
ncdisp(ncname);
%Get file information (and show dimension names and variable names)
finfo=ncinfo(ncname);
dimNames = {finfo.Dimensions.Name};
varNames = {finfo.Variables.Name};
disp(dimNames);
disp(varNames);
%Load the variables
dlon = ncread(ncname,'lon');
dlat = ncread(ncname,'lat');
dtime = ncread(ncname,'time');
CRU_obs_temp = ncread(ncname,'tmp');

% Visually check data by plotting the total period mean
CRU_obs_temp2 = nanmean(CRU_obs_temp,3); 

figure
hold on
axesm miller
worldmap world
geoshow(dlatx,dlonx,CRU_obs_temp2,'displaytype','texturemap');
load coast
colorbar
caxis([-30 30])
title('Annual Mean Temp CRU_obs')


%% Work with MPI Data

% MPI historical
ncname='./../data/tas_Amon_MPI-ESM-LR_historical_r3i1p1_196101-200012_remap.nc';

%Display the metadata for the file (go though to familiarize with the format)
ncdisp(ncname);
%Get file information (and show dimension names and variable names)
finfo=ncinfo(ncname);
dimNames = {finfo.Dimensions.Name};
varNames = {finfo.Variables.Name};
disp(dimNames);
disp(varNames);
%Load the variables
dlon = ncread(ncname,'lon');
dlat = ncread(ncname,'lat');
dtime = ncread(ncname,'time');
MPI_historical_temp = ncread(ncname,'tas');

MPI_historical_temp = MPI_historical_temp - 273.15; % C to K

for i=1:length(dlon)
  dlatx(i,:)=dlat;
end
for i=1:length(dlat)
  dlonx(:,i)=dlon;
end

% Visually check data by plotting the total period mean
MPI_historical_temp2 = nanmean(MPI_historical_temp,3); 

figure
hold on
axesm miller
worldmap world
geoshow(dlatx,dlonx,MPI_historical_temp2,'displaytype','texturemap');
load coast
colorbar
caxis([-30 30])
title('Total period Mean Temperature MPI_historical')


% MPI rcp45
ncname='./../data/tas_Amon_MPI-ESM-LR_rcp45_r3i1p1_207101-210012_remap.nc';

%Display the metadata for the file (go though to familiarize with the format)
ncdisp(ncname);
%Get file information (and show dimension names and variable names)
finfo=ncinfo(ncname);
dimNames = {finfo.Dimensions.Name};
varNames = {finfo.Variables.Name};
disp(dimNames);
disp(varNames);
%Load the variables
dlon = ncread(ncname,'lon');
dlat = ncread(ncname,'lat');
dtime = ncread(ncname,'time');
MPI_rcp45_temp = ncread(ncname,'tas');

MPI_rcp45_temp = MPI_rcp45_temp - 273.15; % C to K

for i=1:length(dlon)
  dlatx(i,:)=dlat;
end
for i=1:length(dlat)
  dlonx(:,i)=dlon;
end

% Visually check data by plotting the total period mean
MPI_rcp45_temp2 = nanmean(MPI_rcp45_temp,3); 

figure
hold on
axesm miller
worldmap world
geoshow(dlatx,dlonx,MPI_rcp45_temp2,'displaytype','texturemap');
load coast
colorbar
caxis([-30 30])
title('Total period Mean Temperature MPI_rcp45')

% MPI rcp85
ncname='./../data/tas_Amon_MPI-ESM-LR_rcp85_r3i1p1_207101-210012_remap.nc';

%Display the metadata for the file (go though to familiarize with the format)
ncdisp(ncname);
%Get file information (and show dimension names and variable names)
finfo=ncinfo(ncname);
dimNames = {finfo.Dimensions.Name};
varNames = {finfo.Variables.Name};
disp(dimNames);
disp(varNames);
%Load the variables
dlon = ncread(ncname,'lon');
dlat = ncread(ncname,'lat');
dtime = ncread(ncname,'time');
MPI_rcp85_temp = ncread(ncname,'tas');

MPI_rcp85_temp = MPI_rcp85_temp - 273.15; % C to K

for i=1:length(dlon)
  dlatx(i,:)=dlat;
end
for i=1:length(dlat)
  dlonx(:,i)=dlon;
end

% Visually check data by plotting the total period mean
MPI_rcp85_temp2 = nanmean(MPI_rcp85_temp,3); 

figure
hold on
axesm miller
worldmap world
geoshow(dlatx,dlonx,MPI_rcp85_temp2,'displaytype','texturemap');
load coast
colorbar
caxis([-30 30])
title('Total period Mean Temperature MPI_rcp45')



%%
% Select subregion (select your own)

% Subregion 1 (total: x: 720, y: 360
First_lon_sub1 = 341; 
Last_lon_sub1 = 435; 
First_lat_sub1 = 250; 
Last_lat_sub1 = 270; 

% Apply to models
CCSM4_historical_temp_subregion1 = CCSM4_historical_temp(First_lon_sub1:Last_lon_sub1,First_lat_sub1:Last_lat_sub1,:);
EC_EARTH_DMI_historical_temp_subregion1 = EC_EARTH_DMI_historical_temp(First_lon_sub1:Last_lon_sub1,First_lat_sub1:Last_lat_sub1,:);
HadGEM2_CC_historical_temp_subregion1 = HadGEM2_CC_historical_temp(First_lon_sub1:Last_lon_sub1,First_lat_sub1:Last_lat_sub1,:);
MPI_historical_temp_subregion1 = MPI_historical_temp(First_lon_sub1:Last_lon_sub1,First_lat_sub1:Last_lat_sub1,:);
%%MPI_rcp45_temp_subregion1 = MPI_rcp45_temp(First_lon_sub1:Last_lon_sub1,First_lat_sub1:Last_lat_sub1,:);
%%MPI_rcp85_temp_subregion1 = MPI_rcp85_temp(First_lon_sub1:Last_lon_sub1,First_lat_sub1:Last_lat_sub1,:);

% Apply to observation data
CRU_obs_temp_subregion1 = CRU_obs_temp(First_lon_sub1:Last_lon_sub1,First_lat_sub1:Last_lat_sub1,:);

% Apply to lat/lon data
dlatx_subregion1 = dlatx(First_lon_sub1:Last_lon_sub1,First_lat_sub1:Last_lat_sub1,:);
dlonx_subregion1 = dlonx(First_lon_sub1:Last_lon_sub1,First_lat_sub1:Last_lat_sub1,:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualize the subregion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CRU_test_subregion_mean = mean(CRU_obs_temp_subregion1,3);


figure
hold on
axesm miller
worldmap world
geoshow(dlatx_subregion1,dlonx_subregion1,CRU_test_subregion_mean,'displaytype','texturemap');
load coast
colorbar
caxis([0 25]) % change to fit region


%%

% This section gives the modelled cells a NAN value in locations where CRU (observed)
% is NAN to be able to meaningfully compare directly

nan_logical = isnan(CRU_obs_temp_subregion1); % check for CRU nan values in selected region

CCSM4_historical_temp_subregion1(nan_logical==1) = nan; % assign nan
EC_EARTH_DMI_historical_temp_subregion1(nan_logical==1) = nan; % assign nan
HadGEM2_CC_historical_temp_subregion1(nan_logical==1) = nan; % assign nan
MPI_historical_temp_subregion1(nan_logical==1) = nan; % assign nan
%MPI_rcp45_temp_subregion1(nan_logical==1) = nan; % assign nan
%MPI_rcp85_temp_subregion1(nan_logical==1) = nan; % assign nan

%%
% % % % % % % % % % % % % % % % % % % % % % % 
% FIGURE 1
% % % % % % % % % % % % % % % % % % % % % % % 

% Take lon mean and squeeze dimension
CCSM4_historical_temp_subregion1_lonmean = nanmean(CCSM4_historical_temp_subregion1,1);
CCSM4_historical_temp_subregion1_lonmean = squeeze(CCSM4_historical_temp_subregion1_lonmean(1,:,:));
EC_EARTH_DMI_historical_temp_subregion1_lonmean = nanmean(EC_EARTH_DMI_historical_temp_subregion1,1);
EC_EARTH_DMI_historical_temp_subregion1_lonmean = squeeze(EC_EARTH_DMI_historical_temp_subregion1_lonmean(1,:,:));
HadGEM2_CC_historical_temp_subregion1_lonmean = nanmean(HadGEM2_CC_historical_temp_subregion1,1);
HadGEM2_CC_historical_temp_subregion1_lonmean = squeeze(HadGEM2_CC_historical_temp_subregion1_lonmean(1,:,:));
CRU_obs_temp_subregion1_lonmean = nanmean(CRU_obs_temp_subregion1,1);
CRU_obs_temp_subregion1_lonmean = squeeze(CRU_obs_temp_subregion1_lonmean(1,:,:));
MPI_historical_temp_subregion1_lonmean = nanmean(MPI_historical_temp_subregion1,1);
%MPI_rcp45_temp_subregion1_lonmean = nanmean(MPI_rcp45_temp_subregion1,1);
%MPI_rcp85_temp_subregion1_lonmean = nanmean(MPI_rcp85_temp_subregion1,1);
% Create weightsnanmean
W = cos(pi/180*dlat); %weights (each cell 0-1)
W_dist = repmat(W,1,720); % distribute weights over all longitudes
% Cut out weights equal to subregion latitudes
W_dist_subregion1 = W_dist(First_lat_sub1:Last_lat_sub1,1);
WSUM_subregion1 = sum(W_dist_subregion1,1); % sum of weights in lat direction for selected region
W_dist_subregion1 = W_dist_subregion1/WSUM_subregion1; % weights to equal 1
% Multiply weights for each of 480 time steps (for 40 years) 

CCSM4_historical_temp_subregion1_lonmean_W = W_dist_subregion1' * CCSM4_historical_temp_subregion1_lonmean;
HadGEM2_CC_historical_temp_subregion1_lonmean_W = W_dist_subregion1'*HadGEM2_CC_historical_temp_subregion1_lonmean;
CRU_obs_temp_subregion1_lonmean_W = W_dist_subregion1' * CRU_obs_temp_subregion1_lonmean;
EC_EARTH_DMI_historical_temp_subregion1_lonmean_W = W_dist_subregion1' * EC_EARTH_DMI_historical_temp_subregion1_lonmean;
%MPI_historical_temp_subregion1_lonmean_W = W_dist_subregion1' * MPI_historical_temp_subregion1_lonmean;
%MPI_rcp45_temp_subregion1_lonmean_W = W_dist_subregion1' * MPI_rcp45_temp_subregion1_lonmean;
%MPI_rcp85_temp_subregion1_lonmean_W = W_dist_subregion1' * MPI_rcp85_temp_subregion1_lonmean;

%TODO:
% CCSM4_historical_temp_subregion1_lonmean_W = times(CCSM4_historical_temp_subregion1_lonmean,W_dist_subregion1');
% EC_EARTH_DMI_historical_temp_subregion1_lonmean_W = EC_EARTH_DMI_historical_temp_subregion1_lonmean.*W_dist_subregion1;
% HadGEM2_CC_historical_temp_subregion1_lonmean_W = HadGEM2_CC_historical_temp_subregion1_lonmean.*W_dist_subregion1;
% CRU_obs_temp_subregion1_lonmean_W = CRU_obs_temp_subregion1_lonmean.*W_dist_subregion1;

% Take sum over the latitudes
% THIS is now the subregion1 weighted 40-year monthly mean time series 
% for the GCM runs
CCSM4_historical_temp_subregion1_lonmean_W = sum(CCSM4_historical_temp_subregion1_lonmean_W,1);
EC_EARTH_DMI_historical_temp_subregion1_lonmean_W = sum(EC_EARTH_DMI_historical_temp_subregion1_lonmean_W,1);
HadGEM2_CC_historical_temp_subregion1_lonmean_W = sum(HadGEM2_CC_historical_temp_subregion1_lonmean_W,1);
CRU_obs_temp_subregion1_lonmean_W = sum(CRU_obs_temp_subregion1_lonmean_W,1);
%MPI_historical_temp_subregion1_lonmean_W = sum(MPI_historical_temp_subregion1_lonmean_W, 1);
%MPI_rcp45_temp_subregion1_lonmean_W = sum(MPI_rcp45_temp_subregion1_lonmean_W, 1);
%MPI_rcp85_temp_subregion1_lonmean_W = sum(MPI_rcp85_temp_subregion1_lonmean_W, 1);

% Set approximate max and min temperatures for your area (see plot)
max_temp = 20;
min_temp = 5;

% employ the quantile-quantile plot to get the X-Y pair data (OBS-SIM)
figure%('visible','off');
hold on
h1 = qqplot(CRU_obs_temp_subregion1_lonmean_W,CCSM4_historical_temp_subregion1_lonmean_W);
h2 = qqplot(CRU_obs_temp_subregion1_lonmean_W,EC_EARTH_DMI_historical_temp_subregion1_lonmean_W);
h3 = qqplot(CRU_obs_temp_subregion1_lonmean_W,HadGEM2_CC_historical_temp_subregion1_lonmean_W);
set(h1(1),'marker','*','markersize',5,'markeredgecolor','b'); % symbols
set(h2(1),'marker','*','markersize',5,'markeredgecolor','g'); % symbols
set(h3(1),'marker','*','markersize',5,'markeredgecolor','r'); % symbols
QQ_data = get(gca,'Children');
X_data = get(QQ_data,'XData'); % Extracts all the X data
Y_data = get(QQ_data,'YData'); % Extracts all the Y data
QQ_CCSM4 = Y_data{7}; % Extracts the specific data from the above 'X_data' og 'Y_data'
QQ_EC_EARTH_DMI = Y_data{4}; % Extracts the specific data from the above 'X_data' og 'Y_data'
QQ_HadGEM2_CC = Y_data{1}; % Extracts the specific data from the above 'X_data' og 'Y_data'
QQ_CRU = X_data{1}; % Extracts the specific data from the above 'X_data' og 'Y_data'
xlim([min_temp max_temp]); % MODIFY according to data but keep the x and y range equal to see the slope
ylim([min_temp max_temp]); % MODIFY according to data but keep the x and y range equal to see the slope
legend
print('-dtiff','-r300','QQ');


% Extract the 50% hottest data
QQ_CRU_50pct = QQ_CRU(241:480);
QQ_CCSM4_50pct = QQ_CCSM4(241:480);
QQ_EC_EARTH_DMI_50pct = QQ_EC_EARTH_DMI(241:480);
QQ_HadGEM2_CC_50pct = QQ_HadGEM2_CC(241:480);

% Perfect fit values
X = [min_temp:0.5:max_temp];
Y = [min_temp:0.5:max_temp];

% X-values to use when plotting trendlines
X_trend = [QQ_CRU_50pct(1):0.25:QQ_CRU_50pct(end)];
% X_trend = [0:0.25:40];

% Slopes of linear regression lines
poly_perfect = polyfit(X,Y,1); 
Slope_PF = poly_perfect(1); % Slope of perfect fit
poly_CCSM4 = polyfit(QQ_CRU_50pct,QQ_CCSM4_50pct,1); 
Slope_CCSM4 = poly_CCSM4(1); % Slope of CCSM4
poly_EC_EARTH_DMI = polyfit(QQ_CRU_50pct,QQ_EC_EARTH_DMI_50pct,1); 
Slope_EC_EARTH_DMI = poly_EC_EARTH_DMI(1); % Slope of EC_EARTH_DMI
poly_HadGEM2_CC = polyfit(QQ_CRU_50pct,QQ_HadGEM2_CC_50pct,1); 
Slope_HadGEM2_CC = poly_HadGEM2_CC(1); % Slope of HadGM2_CC

% Generate the data for plotting the 50pct best fit lines
Y_CCSM4_trendline = poly_CCSM4(1)*X_trend + poly_CCSM4(2);
Y_EC_EARTH_DMI_trendline = poly_EC_EARTH_DMI(1)*X_trend + poly_EC_EARTH_DMI(2);
Y_HadGEM2_CC_trendline = poly_HadGEM2_CC(1)*X_trend + poly_HadGEM2_CC(2);

% The figure below plots the entire data range whereas the slope is only calculated on the basis
% of the hottest 50%
figure
hold on
plot(QQ_CRU,QQ_CCSM4,'*b');
plot(QQ_CRU,QQ_EC_EARTH_DMI,'*g');
plot(QQ_CRU,QQ_HadGEM2_CC,'*r');
plot(X,Y,'-k','LineWidth',1,'Markersize',2,'color',[0.5 0.5 0.5]);
plot(X_trend,Y_CCSM4_trendline,'-k','LineWidth',2);
plot(X_trend,Y_EC_EARTH_DMI_trendline,'-k','LineWidth',2);
plot(X_trend,Y_HadGEM2_CC_trendline,'-k','LineWidth',2);
xlim([min_temp max_temp]);
ylim([min_temp max_temp]);
xlabel('OBS')
ylabel('MODEL')
legend('CCSM4','EC EARTH DMI','HadGEM2 CC','location','northwest');
%print('-dtiff','-r300','PLOT');


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recreate fig 2. From Boberg and Christensen 2012 using the RCP 4.5 data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ncname='./../data/tas_Amon_CCSM4_rcp45_r1i1p1_200601-229912_remap_2071_2100.nc'; % insert correct name

%Display the metadata for the file (go though to familiarize with the format)
ncdisp(ncname);
%Get file information (and show dimension names and variable names)
finfo=ncinfo(ncname);
dimNames = {finfo.Dimensions.Name};
varNames = {finfo.Variables.Name};
disp(dimNames);
disp(varNames);
%Load the variables
dlon = ncread(ncname,'lon');
dlat = ncread(ncname,'lat');
dtime = ncread(ncname,'time');
CCSM4_RCP45_temp = ncread(ncname,'tas'); CCSM4_RCP45_temp = permute(CCSM4_RCP45_temp,[2 1 3]);


CCSM4_RCP45_temp = permute(CCSM4_RCP45_temp,[2 1 3]);
CCSM4_RCP45_temp = flipdim(CCSM4_RCP45_temp,1);
CCSM4_RCP45_temp = CCSM4_RCP45_temp - 273.15; % C to K
% Check that data is rotated correctly by plotting the period mean
CCSM4_RCP45_temp2 = nanmean(CCSM4_RCP45_temp,3); 


figure
hold on
axesm miller
worldmap world
geoshow(dlatx,dlonx,CCSM4_RCP45_temp2,'displaytype','texturemap');
load coast
colorbar
title('OBS: test only to check data orientation and units - mean of TOTAL period'); 
caxis([-30 30])



%%


CCSM4_RCP45_temp_subregion1 = CCSM4_RCP45_temp(First_lon_sub1:Last_lon_sub1,First_lat_sub1:Last_lat_sub1,:); % Select subregion
CCSM4_RCP45_temp_subregion1(nan_logical(:,:,end-359:end)==1) = nan; % assign nan


% Take lon mean and squeeze dimension
CCSM4_RCP45_temp_subregion1_lonmean = nanmean(CCSM4_RCP45_temp_subregion1,1);
CCSM4_RCP45_temp_subregion1_lonmean = squeeze(CCSM4_RCP45_temp_subregion1_lonmean(1,:,:));
% Cut out weights equal to subregion latitudes
W_dist_subregion1 = W_dist(First_lat_sub1:Last_lat_sub1,1);
WSUM_subregion1 = sum(W_dist_subregion1,1); % sum of weights in lat direction for selected region
W_dist_subregion1 = W_dist_subregion1/WSUM_subregion1; % weights to equal 1

% Multiply weights for each of 360 time steps (for 30 years) 
CCSM4_RCP45_temp_subregion1_lonmean_W = CCSM4_RCP45_temp_subregion1_lonmean'*W_dist_subregion1;

% Take sum over the latitudes
CCSM4_RCP45_temp_subregion1_lonmean_W = sum(CCSM4_RCP45_temp_subregion1_lonmean_W,1);

% Sort data and take 50% hottest
CCSM4_RCP45_temp_subregion1_lonmean_W_sort = sort(CCSM4_RCP45_temp_subregion1_lonmean_W);
CCSM4_RCP45_temp_subregion1_lonmean_W_sort_50pct = CCSM4_RCP45_temp_subregion1_lonmean_W_sort(:,181:end);

% take mean of future RCP projection
CCSM4_RCP45_temp_subregion1_mean = nanmean(CCSM4_RCP45_temp_subregion1_lonmean_W_sort_50pct);

% take mean of historial data (previous work)
CCSM4_historical_temp_subregion1_mean = nanmean(QQ_CCSM4_50pct);

% Temp_change (figure y-axis)
CCSM4_temp_change = CCSM4_RCP45_temp_subregion1_mean - CCSM4_historical_temp_subregion1_mean;

% Temp bias slope (x-axis) - all for 50% hottest timesteps
% I.e.: 
% model_slope (T_obs/T_model) / T_bias / T_obs
% Or:
% model_slope (T_obs/T_model) / (T_model - T_obs) / T_obs

CCSM4_temp_bias_slope = Slope_CCSM4/(mean(QQ_CCSM4_50pct)-mean(QQ_CRU_50pct))/mean(QQ_CRU_50pct);

figure
hold on
plot(CCSM4_temp_bias_slope,CCSM4_temp_change,'*k')
% insert data from more models her...
xlabel('bias slope (model_slope/T_bias/T_obs) for the 50% warmest months')
ylabel('Projected temperature change for the 50% warmest months')
title('Fig2, Christensen and Boberg (2012), RCP4.5 2071-2100')









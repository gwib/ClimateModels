%%
run('loadData')
run('conversion')
%% stuff copied from other exercises

%% sum of "anomalies" over region
% EXERCISE 3-4 (W1)

clearvars *lon_mean*

% % CRU

% Code to geometrically weigh the latitudinal direction.

W = cos(pi/180*dlat_CRU)'; %weights (each cell 0-1)
WSUM = sum(W); % sum of weights
W = W/WSUM; % weights to equal 1
sum(W) % Check that they equal 1

no_timesteps_CRU = numel(TEMP_anomaly_CRU(1,1,:));

for i = 1:no_timesteps_CRU
lon_mean_CRU(:,:,i) = nanmean(TEMP_anomaly_CRU(:,:,i)); % lon mean
end
lon_mean_CRU = squeeze(lon_mean_CRU(1,:,:)); % remove dimension

W_time = repmat(W,no_timesteps_CRU,1)'; % repeat weigths over time

lon_mean_CRU_W = W_time.*lon_mean_CRU; % multiply weights with lon means
lon_mean_CRU_W_sum = nansum(lon_mean_CRU_W); % Summarize over latitudes
time_CRU = [1850:1/12:2014+2/12]; % create simple time array for plot


% % GISS

W = cos(pi/180*dlat_GISS)'; %weights (each cell 0-1)
WSUM = sum(W); % sum of weights
W = W/WSUM; % weights to equal 1
sum(W) % Check that they equal 1

no_timesteps_GISS = numel(TEMP_anomaly_GISS(1,1,:));

for i = 1:no_timesteps_GISS
lon_mean_GISS(:,:,i) = nanmean(TEMP_anomaly_GISS(:,:,i)); % lon mean
end
lon_mean_GISS = squeeze(lon_mean_GISS(1,:,:)); % remove dimension

W_time = repmat(W,no_timesteps_GISS,1)'; % repeat weigths over time

lon_mean_GISS_W = W_time.*lon_mean_GISS; % multiply weights with lon means
lon_mean_GISS_W_sum = nansum(lon_mean_GISS_W); % Summarize over latitudes
time_GISS = [1880:1/12:2017+2/12]; % create simple time array for plot




figure
hold on
plot(time_CRU,lon_mean_CRU_W_sum);
plot(time_GISS,lon_mean_GISS_W_sum);
title('CRU and GISS near-surface temperature anomalies');
ylabel('temp (anomaly) (C^o)');
xlabel('year');
h_legend = legend('CRU','GISS','Location','NorthWest'); % Legend text and location
set(h_legend,'Fontsize',11,'Fontweight','bold'); % Legend characteristics
print('-dtiff','-r300','CRU_GISS_temp_anomaly'); % save WITH legend
%%

% worldmap -region-


%% take out values for certain months (W1)


% % % % % % % % 
% % GISS
% % % % % % % % 

clearvars *JUL* *JAN* -except *CRU_JULY* *CRU_JAN* *GISS_JULY* *GISS_JAN*


JANUARY_1 = ((1895-1880)*12)+1;
for i = 1:40
TEMP_anomaly_JAN1(:,:,i) = TEMP_anomaly_GISS(:,:,JANUARY_1);
JANUARY_1 = JANUARY_1+12;
end
TEMP_anomaly_GISS_JAN1 = nanmean(TEMP_anomaly_JAN1,3)';

JANUARY_2 = ((1935-1880)*12)+1;
for i = 1:40
TEMP_anomaly_JAN2(:,:,i) = TEMP_anomaly_GISS(:,:,JANUARY_2);
JANUARY_2 = JANUARY_2+12;
end
TEMP_anomaly_GISS_JAN2 = nanmean(TEMP_anomaly_JAN2,3)';

JANUARY_3 = ((1975-1880)*12)+1;
for i = 1:40
TEMP_anomaly_JAN3(:,:,i) = TEMP_anomaly_GISS(:,:,JANUARY_3);
JANUARY_3 = JANUARY_3+12;
end
TEMP_anomaly_GISS_JAN3 = nanmean(TEMP_anomaly_JAN3,3)';





JULY_1 = ((1895-1880)*12)+7;
for i = 1:40
TEMP_anomaly_JULY1(:,:,i) = TEMP_anomaly_GISS(:,:,JULY_1);
JULY_1 = JULY_1+12;
end
TEMP_anomaly_GISS_JULY1 = nanmean(TEMP_anomaly_JULY1,3)';

JULY_2 = ((1935-1880)*12)+7;
for i = 1:40
TEMP_anomaly_JULY2(:,:,i) = TEMP_anomaly_GISS(:,:,JULY_2);
JULY_2 = JULY_2+12;
end
TEMP_anomaly_GISS_JULY2 = nanmean(TEMP_anomaly_JULY2,3)';

JULY_3 = ((1975-1880)*12)+7;
for i = 1:40
TEMP_anomaly_JULY3(:,:,i) = TEMP_anomaly_GISS(:,:,JULY_3);
JULY_3 = JULY_3+12;
end
TEMP_anomaly_GISS_JULY3 = nanmean(TEMP_anomaly_JULY3,3)';




%To plot a 2D matrix with geoshow the latitudes and longitudes must be converted to matrixes
for i=1:length(dlat_CRU)
  dlonx_GISS(i,:)=dlon_CRU;
end

for i=1:length(dlon_CRU)
  dlatx_GISS(:,i)=dlat_CRU;
end


clearvars *JUL* *JAN* -except *CRU_JULY* *CRU_JAN* *GISS_JULY* *GISS_JAN*


%%
% total period mean (W3)
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



%% Plot for subregion (W3)
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




% This section gives the modelled cells a NAN value in locations where CRU (observed)
% is NAN to be able to meaningfully compare directly

nan_logical = isnan(CRU_obs_temp_subregion1); % check for CRU nan values in selected region

CCSM4_historical_temp_subregion1(nan_logical==1) = nan; % assign nan
EC_EARTH_DMI_historical_temp_subregion1(nan_logical==1) = nan; % assign nan
HadGEM2_CC_historical_temp_subregion1(nan_logical==1) = nan; % assign nan
MPI_historical_temp_subregion1(nan_logical==1) = nan; % assign nan
%MPI_rcp45_temp_subregion1(nan_logical==1) = nan; % assign nan
%MPI_rcp85_temp_subregion1(nan_logical==1) = nan; % assign nan


% % % % % % % % % % % % % % % % % % % % % % % 
% FIGURE 1 -- correlation plots
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

%% Map comparison between datasets (W2)

%first we read the files

%Load the variables
dlon = ncread(ncname,'longitude');
dlat = ncread(ncname,'latitude');
dtime = ncread(ncname,'time');
EOBS_TEMP = ncread(ncname,'tg');

% Generate distributed lats and lons
for i=1:length(dlon)
  dlatx(i,:)=dlat;
end
for i=1:length(dlat)
  dlonx(:,i)=dlon;
end

% Check that data is rotated correctly by plotting the period mean
Check_EOBS_TEMP = nanmean(EOBS_TEMP,3); 

figure
hold on
worldmap europe
colorbar
geoshow(dlatx,dlonx,Check_EOBS_TEMP,'displaytype','texturemap');
load coastlines
geoshow(coastlat,coastlon,'color','black')


%final plots
figure

subplot(2,2,1)
hold on
axesm miller
worldmap europe
geoshow(dlatx,dlonx,Check_EOBS_TEMP,'displaytype','texturemap');
load coast
colorbar
caxis([-5 25])
title('EOBS TEMP')
load coastlines
geoshow(coastlat,coastlon,'color','black')

subplot(2,2,2)
hold on
axesm miller
worldmap europe
geoshow(dlatx,dlonx,Check_ERAI_TEMP,'displaytype','texturemap');
load coast
colorbar
caxis([-5 25])
title('ERAI TEMP')
load coastlines
geoshow(coastlat,coastlon,'color','black')

subplot(2,2,3)
hold on
axesm miller
worldmap europe
geoshow(dlatx,dlonx,Check_NCEP_NCAR_TEMP,'displaytype','texturemap');
load coast
colorbar
caxis([-5 25])
title('NCEP TEMP')
load coastlines
geoshow(coastlat,coastlon,'color','black')


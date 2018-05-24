ncname='./../data/CERES_SWUPALL_200003-201012.nc'; % insert correct name

%Display the metadata for the file (go though to familiarize with the format)
ncdisp(ncname);

%Get file information (and show dimension names and variable names)
finfo=ncinfo(ncname);
dimNames = {finfo.Dimensions.Name};
varNames = {finfo.Variables.Name};
disp(dimNames);
disp(varNames);

%Load the variables
dlon=ncread(ncname,'lon');
dlat=ncread(ncname,'lat');
dtime=ncread(ncname,'time');
dSW=ncread(ncname,'toa_sw_all_mon');

%%
% Test to see data

%Plot the timeseries for a single grid point

lon_cell = 50;
lat_cell = 50;

dSWt(1:130)=dSW(lon_cell,lat_cell,1:130); % extract time period
figure; 
plot(dtime,dSWt);

%%

% Assignment 1

JAN2005 = 71;
JUL2005 = 77;

%Extract data from January 2005 
dSW_JAN(:,:)=dSW(:,:,JAN2005);
%and transpose the matrix
dSW_JAN=dSW_JAN';

%Extract data from July 2005 
dSW_JUL(:,:)=dSW(:,:,JUL2005);
%and transpose the matrix
dSW_JUL=dSW_JUL';

%To plot a 2D matrix with geoshow the latitudes and longitudes must be converted to matrixes
for i=1:length(dlat)
  dlonx(i,:)=dlon;
end

for i=1:length(dlon)
  dlatx(:,i)=dlat;
end

%Plot the data;
% JANUARY
figure%('visible','off');
worldmap world
geoshow(dlatx,dlonx,dSW_JAN,'displaytype','texturemap');
% contourcmap('Colorbar', 'on','Location', 'horizontal');
load coast
geoshow(lat,long,'color','w')
colorbar;
caxis([0 325])
title('SW up for January 2005')
print('-dtiff','-r450','SW_JAN');
% JULY
figure%('visible','off');
worldmap world
geoshow(dlatx,dlonx,dSW_JUL,'displaytype','texturemap');
% contourcmap('Colorbar', 'on','Location', 'horizontal');
load coast
geoshow(lat,long,'color','w')
colorbar;
caxis([0 325])
title('SW up for July 2005')
print('-dtiff','-r450','SW_JUL');

%%

% Assignment 2
W = cos(pi/180*dlat)'; %weights (each cell 0-1)
WSUM = sum(W); % sum of weights
W = W/WSUM; % weights to equal 1
sum(W)

for i = 1:130
lon_mean(:,:,i) = mean(dSW(:,:,i)); % lon mean
end
lon_mean = squeeze(lon_mean(1,:,:)); % remove dimension

W_time = repmat(W,130,1)'; % repeat weigths over time

lon_mean_W = W_time.*lon_mean; % multiply weights with lon means
lon_mean_W_sum = sum(lon_mean_W); % Summarize over latitudes
time = [1:1:130]; % create simple time array for plot


figure
hold on
plot(time,lon_mean_W_sum);

%%

% Assignment 3

SW_UP_2001_2010_mean = mean(lon_mean_W_sum(11:130));
SW_down = 1366/4;
Net_SW = SW_down - SW_UP_2001_2010_mean;

%%

% % Assignment 4

% From exercise 1 mean(OLR_CERES_lon_mean_W_SUM)
MEAN_OLR = 239.8109;

TOA_CERES_energy_balance = Net_SW - MEAN_OLR;

fprintf('CERES energy balance: %.4f', TOA_CERES_energy_balance)














% 
% 
% 
% 
% %Create a latitudinal mean (zonal mean) by making a mean row wise (second dimension)
% %Here for the first day:
% clear dolr1
% dolr1(:,:)=dolr(:,:,1);
% dolr1=dolr1';
% meanb=mean(dolr1,2);
% 
% figure; 
% plot(dlat,meanb);
% 
% %A time average is done by taking the mean over the time dimension (3rd dimension)
% meant=mean(dolr,3);
% meant=meant';
% 
% figure;
% worldmap world
% geoshow(dlatx,dlonx,meant,'displaytype','texturemap');
% colorbar;
% load coast
% geoshow(lat,long,'color','w')

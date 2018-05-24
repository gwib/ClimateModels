ncname='./../data/gistemp1200_ERSSTv4.nc'; % insert correct name

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
dtime_bnds = ncread(ncname,'time_bnds');
dtempanomaly = ncread(ncname,'tempanomaly');


% start day: 15/01/1880 (1/1/1880 + min(dtime))
% end day: 15/03/2017 (1/1/1880 + max(dtime))
%%

First_30 = dtempanomaly(:,:,1:360); % Earliest period
Last_30 = dtempanomaly(:,:,end-359:end); % Latest period

First_30_mean = nanmean(First_30,3);
Last_30_mean = nanmean(Last_30,3);

Diff = (Last_30_mean - First_30_mean)';

%To plot a 2D matrix with geoshow the latitudes and longitudes must be converted to matrixes
for i=1:length(dlat)
  dlonx(i,:)=dlon;
end

for i=1:length(dlon)
  dlatx(:,i)=dlat;
end

figure%('visible','off'); % Open the figure from your work directory!
worldmap world
geoshow(dlatx,dlonx,Diff,'displaytype','texturemap');
% contourcmap('Colorbar', 'on','Location', 'horizontal');
load coast
geoshow(lat,long,'color','w')
colorbar;
newmap = jet;                    %starting map
ncol = size(newmap,1);           %how big is it?
zpos = 1 + floor(0.01 * ncol);    % 0.01 of way through
newmap(zpos,:) = [1 1 1];        %set that position to white
colormap(newmap);                %activate it
print('-dtiff','-r300','Anomaly');


 
%%

% Make plot of all years

All_lonmean = nanmean(dtempanomaly,1); % lon mean all
All_lonmean = squeeze(All_lonmean(1,:,:)); % remove dim latest

No_years = round(numel(dtempanomaly(1,1,:))/12); % please understand this line

% Make yearly means total period
k = 1;
for i = 1:No_years
    All_YEAR(:,i) = nanmean(All_lonmean(:,k:k+11),2);
    k = k+12;
end

W = cos(pi/180*dlat)'; %weights (each cell 0-1)
WSUM = sum(W); % sum of weights
W = W/WSUM; % normalize weights
SUM_check = sum(W); % Check that the weights equal 1
W_time_YEAR_30 = repmat(W,30,1)'; % repeat weigths over time - yearly
W_time_YEAR_ALL = repmat(W,No_years,1)'; % repeat weigths over time - yearly
All_W_YEAR = W_time_YEAR_ALL.*All_YEAR; % Apply weights
All_W_YEAR_latsum = nansum(All_W_YEAR); % Summarize over latitudes
time_YEAR_ALL = [1:1:No_years]; % create simple time array for plot - YEAR

% % Make plot

x_axis = [1880:1:2016];

figure
hold on
plot(x_axis,All_W_YEAR_latsum);
xlabel('year')
ylabel('anomaly (deg K/C)')
title('Entire period anomaly');
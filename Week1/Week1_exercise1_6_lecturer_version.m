clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WEEK1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% EXERCISE1

% CRU

ncname_CRU='./../data/CRUTEM.4.2.0.0.anomalies.nc';
%Display the metadata for the file (go though to familiarize with the format)
ncdisp(ncname_CRU);

%Get file information (and show dimension names and variable names)
finfo=ncinfo(ncname_CRU);
dimNames = {finfo.Dimensions.Name};
varNames = {finfo.Variables.Name};
disp(dimNames);
disp(varNames);

%Load the variables
dlon_CRU = ncread(ncname_CRU,'longitude');
dlat_CRU = ncread(ncname_CRU,'latitude');
dtime_CRU = ncread(ncname_CRU,'time');
TEMP_anomaly_CRU = ncread(ncname_CRU,'temperature_anomaly');

% GISS

ncname_GISS='./../data/gistemp1200_ERSSTv4.nc';
%Display the metadata for the file (go though to familiarize with the format)
ncdisp(ncname_GISS);

%Get file ncname_GISS (and show dimension names and variable names)
finfo=ncinfo(ncname_GISS);
dimNames = {finfo.Dimensions.Name};
varNames = {finfo.Variables.Name};
disp(dimNames);
disp(varNames);

%Load the variables
dlon_GISS = ncread(ncname_GISS,'lon');
dlat_GISS = ncread(ncname_GISS,'lat');
dtime_GISS = ncread(ncname_GISS,'time');
TEMP_anomaly_GISS = ncread(ncname_GISS,'tempanomaly');

%%

% EXERCISE2

% plot mean JANUARY anomaly values for each of the four latest 
% 40-year time slices:
% 1855-1894 - 1895-1934 - 1935-1974  1975-2014

% % % % % % % 
% % % % % CRU
% % % % % % % 

clearvars JUL* JAN*

JANUARY_1 = ((1855-1850)*12)+1; % why are we doing this??

%TEMP_anomaly_JAN1 = zeros(36,72); % preallocation for performance
%disp(TEMP_anomaly_JAN1)
for i = 1:30
TEMP_anomaly_JAN1(:,:,i) = TEMP_anomaly_CRU(:,:,JANUARY_1);
JANUARY_1 = JANUARY_1+12;
end
TEMP_anomaly_CRU_JAN1 = nanmean(TEMP_anomaly_JAN1,3)';
dispc(TEMP_anomaly_CRU_JAN1)
JANUARY_2 = ((1895-1850)*12)+1;
TEMP_anomaly_JAN2 = zeros(36,72); % preallocation for performance

for i = 1:40
TEMP_anomaly_JAN2(:,:,i) = TEMP_anomaly_CRU(:,:,JANUARY_2);
JANUARY_2 = JANUARY_2+12;
end
TEMP_anomaly_CRU_JAN2 = nanmean(TEMP_anomaly_JAN2,3)';

JANUARY_3 = ((1935-1850)*12)+1;
TEMP_anomaly_JAN3 = zeros(36,72); % preallocation for performance

for i = 1:40
TEMP_anomaly_JAN3(:,:,i) = TEMP_anomaly_CRU(:,:,JANUARY_3);
JANUARY_3 = JANUARY_3+12;
end
TEMP_anomaly_CRU_JAN3 = nanmean(TEMP_anomaly_JAN3,3)';

JANUARY_4 = ((1975-1850)*12)+1;
TEMP_anomaly_JAN4 = zeros(36,72); % preallocation for performance

for i = 1:40
TEMP_anomaly_JAN4(:,:,i) = TEMP_anomaly_CRU(:,:,JANUARY_4);
JANUARY_4 = JANUARY_4+12;
end
TEMP_anomaly_CRU_JAN4 = nanmean(TEMP_anomaly_JAN4,3)'; % Notice that the ' transposes the data






JULY_1 = ((1855-1850)*12)+7;
TEMP_anomaly_JULY1 = zeros(36,72); % preallocation for performance

for i = 1:40
TEMP_anomaly_JULY1(:,:,i) = TEMP_anomaly_CRU(:,:,JULY_1);
JULY_1 = JULY_1+12;
end
TEMP_anomaly_CRU_JULY1 = nanmean(TEMP_anomaly_JULY1,3)';

JULY_2 = ((1895-1850)*12)+7;
TEMP_anomaly_JULY2 = zeros(36,72); % preallocation for performance

for i = 1:40
TEMP_anomaly_JULY2(:,:,i) = TEMP_anomaly_CRU(:,:,JULY_2);
JULY_2 = JULY_2+12;
end
TEMP_anomaly_CRU_JULY2 = nanmean(TEMP_anomaly_JULY2,3)';

JULY_3 = ((1935-1850)*12)+7;
TEMP_anomaly_JULY3 = zeros(36,72); % preallocation for performance

for i = 1:40
TEMP_anomaly_JULY3(:,:,i) = TEMP_anomaly_CRU(:,:,JULY_3);
JULY_3 = JULY_3+12;
end
TEMP_anomaly_CRU_JULY3 = nanmean(TEMP_anomaly_JULY3,3)';

JULY_4 = ((1975-1850)*12)+7;
TEMP_anomaly_JULY4 = zeros(36,72); % preallocation for performance

for i = 1:39
TEMP_anomaly_JULY4(:,:,i) = TEMP_anomaly_CRU(:,:,JULY_4);
JULY_4 = JULY_4+12;
end
TEMP_anomaly_CRU_JULY4 = nanmean(TEMP_anomaly_JULY4,3)'; % Notice that the ' transposes the data



%To plot a 2D matrix with geoshow the latitudes and longitudes must be converted to matrixes
for i=1:length(dlat_CRU)
  dlonx_CRU(i,:)=dlon_CRU;
end

for i=1:length(dlon_CRU)
  dlatx_CRU(:,i)=dlat_CRU;
end


%%


% % % % % % % % 
% % GISS
% % % % % % % % 

clearvars *JUL* *JAN* -except *CRU_JULY* *CRU_JAN* *GISS_JULY* *GISS_JAN*


JANUARY_1 = ((1895-1880)*12)+1;
for i = 1:40
TEMP_anomaly_JAN1(:,:,i) = TEMP_anomaly_GISS(:,:,JANUARY_1);
JANUARY_1 = JANUARY_1+12
end
TEMP_anomaly_GISS_JAN1 = nanmean(TEMP_anomaly_JAN1,3)';

JANUARY_2 = ((1935-1880)*12)+1
for i = 1:40
TEMP_anomaly_JAN2(:,:,i) = TEMP_anomaly_GISS(:,:,JANUARY_2);
JANUARY_2 = JANUARY_2+12
end
TEMP_anomaly_GISS_JAN2 = nanmean(TEMP_anomaly_JAN2,3)';

JANUARY_3 = ((1975-1880)*12)+1
for i = 1:40
TEMP_anomaly_JAN3(:,:,i) = TEMP_anomaly_GISS(:,:,JANUARY_3);
JANUARY_3 = JANUARY_3+12
end
TEMP_anomaly_GISS_JAN3 = nanmean(TEMP_anomaly_JAN3,3)';





JULY_1 = ((1895-1880)*12)+7
for i = 1:40
TEMP_anomaly_JULY1(:,:,i) = TEMP_anomaly_GISS(:,:,JULY_1);
JULY_1 = JULY_1+12
end
TEMP_anomaly_GISS_JULY1 = nanmean(TEMP_anomaly_JULY1,3)';

JULY_2 = ((1935-1880)*12)+7
for i = 1:40
TEMP_anomaly_JULY2(:,:,i) = TEMP_anomaly_GISS(:,:,JULY_2);
JULY_2 = JULY_2+12
end
TEMP_anomaly_GISS_JULY2 = nanmean(TEMP_anomaly_JULY2,3)';

JULY_3 = ((1975-1880)*12)+7
for i = 1:40
TEMP_anomaly_JULY3(:,:,i) = TEMP_anomaly_GISS(:,:,JULY_3);
JULY_3 = JULY_3+12
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

% % % % % CRU

% JAN

figure('visible','off');
hold on

subplot(2,2,1)
set(gca,'Position',[0 0.5 0.5 0.5])
worldmap world
% worldmap([-50 50],[160 -30])
geoshow(dlatx_CRU,dlonx_CRU,TEMP_anomaly_CRU_JAN1,'displaytype','texturemap');
% contourcmap('Colorbar', 'on','Location', 'horizontal');
load coast
geoshow(lat,long,'color','w'); colorbar;
caxis([-3 3])
% Make sure NAN-values are white:
newmap = jet;                    %starting map
ncol = size(newmap,1);           %how big is it?
zpos = 1 + floor(0.01 * ncol);    % 0.01 of way through
newmap(zpos,:) = [1 1 1];        %set that position to white
colormap(newmap);                %activate it
title('JAN 1855-1894');

subplot(2,2,2)
set(gca,'Position',[0.5 0.5 0.5 0.5])
worldmap world
% worldmap europe
geoshow(dlatx_CRU,dlonx_CRU,TEMP_anomaly_CRU_JAN2,'displaytype','texturemap');
% contourcmap('Colorbar', 'on','Location', 'horizontal');
load coast
geoshow(lat,long,'color','w'); colorbar;
caxis([-3 3])
% Make sure NAN-values are white:
newmap = jet;                    %starting map
ncol = size(newmap,1);           %how big is it?
zpos = 1 + floor(0.01 * ncol);    % 0.01 of way through
newmap(zpos,:) = [1 1 1];        %set that position to white
colormap(newmap);                %activate it
title('JAN 1895-1934');

subplot(2,2,3)
set(gca,'Position',[0 0 0.5 0.5])
worldmap world
geoshow(dlatx_CRU,dlonx_CRU,TEMP_anomaly_CRU_JAN3,'displaytype','texturemap');
% contourcmap('Colorbar', 'on','Location', 'horizontal');
load coast
geoshow(lat,long,'color','w'); colorbar;
caxis([-3 3])
% Make sure NAN-values are white:
newmap = jet;                    %starting map
ncol = size(newmap,1);           %how big is it?
zpos = 1 + floor(0.01 * ncol);    % 0.01 of way through
newmap(zpos,:) = [1 1 1];        %set that position to white
colormap(newmap);                %activate it
title('JAN 1935-1974');

subplot(2,2,4)
set(gca,'Position',[0.5 0 0.5 0.5])
worldmap world
geoshow(dlatx_CRU,dlonx_CRU,TEMP_anomaly_CRU_JAN4,'displaytype','texturemap');
% contourcmap('Colorbar', 'on','Location', 'horizontal');
load coast
geoshow(lat,long,'color','w'); colorbar;
caxis([-3 3])
% Make sure NAN-values are white:
newmap = jet;                    %starting map
ncol = size(newmap,1);           %how big is it?
zpos = 1 + floor(0.01 * ncol);    % 0.01 of way through
newmap(zpos,:) = [1 1 1];        %set that position to white
colormap(newmap); %activate it
title('JAN 1975-2014');

print('-dtiff','-r300','TEMP_anomaly_CRU_JAN_fig');


% JUL

figure('visible','off');
hold on

subplot(2,2,1)
set(gca,'Position',[0 0.5 0.5 0.5])
worldmap world
% worldmap([-50 50],[160 -30])
geoshow(dlatx_CRU,dlonx_CRU,TEMP_anomaly_CRU_JULY1,'displaytype','texturemap');
% contourcmap('Colorbar', 'on','Location', 'horizontal');
load coast
geoshow(lat,long,'color','w'); colorbar;
caxis([-3 3])
% Make sure NAN-values are white:
newmap = jet;                    %starting map
ncol = size(newmap,1);           %how big is it?
zpos = 1 + floor(0.01 * ncol);    % 0.01 of way through
newmap(zpos,:) = [1 1 1];        %set that position to white
colormap(newmap);                %activate it
title('JULY 1855-1894');

subplot(2,2,2)
set(gca,'Position',[0.5 0.5 0.5 0.5])
worldmap world
% worldmap europe
geoshow(dlatx_CRU,dlonx_CRU,TEMP_anomaly_CRU_JULY2,'displaytype','texturemap');
% contourcmap('Colorbar', 'on','Location', 'horizontal');
load coast
geoshow(lat,long,'color','w'); colorbar;
caxis([-3 3])
% Make sure NAN-values are white:
newmap = jet;                    %starting map
ncol = size(newmap,1);           %how big is it?
zpos = 1 + floor(0.01 * ncol);    % 0.01 of way through
newmap(zpos,:) = [1 1 1];        %set that position to white
colormap(newmap);                %activate it
title('JULY 1895-1934');

subplot(2,2,3)
set(gca,'Position',[0 0 0.5 0.5])
worldmap world
geoshow(dlatx_CRU,dlonx_CRU,TEMP_anomaly_CRU_JULY3,'displaytype','texturemap');
% contourcmap('Colorbar', 'on','Location', 'horizontal');
load coast
geoshow(lat,long,'color','w'); colorbar;
caxis([-3 3])
% Make sure NAN-values are white:
newmap = jet;                    %starting map
ncol = size(newmap,1);           %how big is it?
zpos = 1 + floor(0.01 * ncol);    % 0.01 of way through
newmap(zpos,:) = [1 1 1];        %set that position to white
colormap(newmap);                %activate it
title('JULY 1935-1974');

subplot(2,2,4)
set(gca,'Position',[0.5 0 0.5 0.5])
worldmap world
geoshow(dlatx_CRU,dlonx_CRU,TEMP_anomaly_CRU_JULY4,'displaytype','texturemap');
% contourcmap('Colorbar', 'on','Location', 'horizontal');
load coast
geoshow(lat,long,'color','w'); colorbar;
caxis([-3 3])
% Make sure NAN-values are white:
newmap = jet;                    %starting map
ncol = size(newmap,1);           %how big is it?
zpos = 1 + floor(0.01 * ncol);    % 0.01 of way through
newmap(zpos,:) = [1 1 1];        %set that position to white
colormap(newmap); %activate it
title('JULY 1975-2013 (39 years)');

print('-dtiff','-r300','TEMP_anomaly_CRU_JULY_fig');




% % GISS

% JAN


figure('visible','off');
hold on

subplot(2,2,1)


subplot(2,2,2)
set(gca,'Position',[0.5 0.5 0.5 0.5])
worldmap world
% worldmap europe
geoshow(dlatx_GISS,dlonx_GISS,TEMP_anomaly_GISS_JAN1,'displaytype','texturemap');
% contourcmap('Colorbar', 'on','Location', 'horizontal');
load coast
geoshow(lat,long,'color','w'); colorbar;
caxis([-3 3])
% Make sure NAN-values are white:
newmap = jet;                    %starting map
ncol = size(newmap,1);           %how big is it?
zpos = 1 + floor(0.01 * ncol);    % 0.01 of way through
newmap(zpos,:) = [1 1 1];        %set that position to white
colormap(newmap);                %activate it
title('JAN 1895-1934');

subplot(2,2,3)
set(gca,'Position',[0 0 0.5 0.5])
worldmap world
geoshow(dlatx_GISS,dlonx_GISS,TEMP_anomaly_GISS_JAN2,'displaytype','texturemap');
% contourcmap('Colorbar', 'on','Location', 'horizontal');
load coast
geoshow(lat,long,'color','w'); colorbar;
caxis([-3 3])
% Make sure NAN-values are white:
newmap = jet;                    %starting map
ncol = size(newmap,1);           %how big is it?
zpos = 1 + floor(0.01 * ncol);    % 0.01 of way through
newmap(zpos,:) = [1 1 1];        %set that position to white
colormap(newmap);                %activate it
title('JAN 1935-1974');

subplot(2,2,4)
set(gca,'Position',[0.5 0 0.5 0.5])
worldmap world
geoshow(dlatx_GISS,dlonx_GISS,TEMP_anomaly_GISS_JAN3,'displaytype','texturemap');
% contourcmap('Colorbar', 'on','Location', 'horizontal');
load coast
geoshow(lat,long,'color','w'); colorbar;
caxis([-3 3])
% Make sure NAN-values are white:
newmap = jet;                    %starting map
ncol = size(newmap,1);           %how big is it?
zpos = 1 + floor(0.01 * ncol);    % 0.01 of way through
newmap(zpos,:) = [1 1 1];        %set that position to white
colormap(newmap); %activate it
title('JAN 1975-2014');

print('-dtiff','-r300','TEMP_anomaly_GISS_JAN_fig');


% JUL

figure('visible','off');
hold on

subplot(2,2,1)


subplot(2,2,2)
set(gca,'Position',[0.5 0.5 0.5 0.5])
worldmap world
% worldmap europe
geoshow(dlatx_GISS,dlonx_GISS,TEMP_anomaly_GISS_JULY1,'displaytype','texturemap');
% contourcmap('Colorbar', 'on','Location', 'horizontal');
load coast
geoshow(lat,long,'color','w'); colorbar;
caxis([-3 3])
% Make sure NAN-values are white:
newmap = jet;                    %starting map
ncol = size(newmap,1);           %how big is it?
zpos = 1 + floor(0.01 * ncol);    % 0.01 of way through
newmap(zpos,:) = [1 1 1];        %set that position to white
colormap(newmap);                %activate it
title('JULY 1895-1934');

subplot(2,2,3)
set(gca,'Position',[0 0 0.5 0.5])
worldmap world
geoshow(dlatx_GISS,dlonx_GISS,TEMP_anomaly_GISS_JULY2,'displaytype','texturemap');
% contourcmap('Colorbar', 'on','Location', 'horizontal');
load coast
geoshow(lat,long,'color','w'); colorbar;
caxis([-3 3])
% Make sure NAN-values are white:
newmap = jet;                    %starting map
ncol = size(newmap,1);           %how big is it?
zpos = 1 + floor(0.01 * ncol);    % 0.01 of way through
newmap(zpos,:) = [1 1 1];        %set that position to white
colormap(newmap);                %activate it
title('JULY 1935-1974');

subplot(2,2,4)
set(gca,'Position',[0.5 0 0.5 0.5])
worldmap world
geoshow(dlatx_GISS,dlonx_GISS,TEMP_anomaly_GISS_JULY3,'displaytype','texturemap');
% contourcmap('Colorbar', 'on','Location', 'horizontal');
load coast
geoshow(lat,long,'color','w'); colorbar;
caxis([-3 3])
% Make sure NAN-values are white:
newmap = jet;                    %starting map
ncol = size(newmap,1);           %how big is it?
zpos = 1 + floor(0.01 * ncol);    % 0.01 of way through
newmap(zpos,:) = [1 1 1];        %set that position to white
colormap(newmap); %activate it
title('JULY 1975-2014');

print('-dtiff','-r300','TEMP_anomaly_GISS_JULY_fig');



%%

% EXERCISE 3-4

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


% EXERCISE 5...




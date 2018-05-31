clear all

latLim = [58, 84];
lonLim = [-91, 29];

ncname1='./data/absolute_remap.nc'; % insert correct name
ncname2='./data/HadCRUT.4.6.0.0.median_remap.nc'; % insert correct name
ncdisp(ncname1);
finfo=ncinfo(ncname1);
dimNames = {finfo.Dimensions.Name};
varNames = {finfo.Variables.Name};
disp(dimNames);
disp(varNames);
abs = ncread(ncname1,'tem');
lat = ncread(ncname1,'lat');
lon = ncread(ncname1,'lon');
dtime_hadCrut = ncread(ncname2, 'time');
%To plot a 2D matrix with geoshow the latitudes and longitudes must be converted to matrixes
for i=1:length(lat)
  dlonx(i,:)=lon;
end
dlonx = dlonx';
for i=1:length(lon)
  dlatx(:,i)=lat;
end
dlatx = dlatx';
ncdisp(ncname2);
finfo=ncinfo(ncname2);
dimNames = {finfo.Dimensions.Name};
varNames = {finfo.Variables.Name};
disp(dimNames);
disp(varNames);
anom = ncread(ncname2,'temperature_anomaly');
% m = repmat([1:12],1,168); m = cat(2,m,[1 2 3]);
abs2 = repmat(abs,1,1,168); %168 full years
abs2 = cat(3,abs2,abs(:,:,1),abs(:,:,2),abs(:,:,3));
for i = 1:2019
    anom2(:,:,i) = anom(:,:,i)+abs2(:,:,i);
end
anom2 = anom2(:,:,1561:1872);
anom2_mean = nanmean(anom2,3);

 
 
figure
hold on
worldmap world
geoshow(dlatx,dlonx,anom2_mean,'displaytype','texturemap');
colorbar
caxis([-30 30])
load coastlines
geoshow(coastlat,coastlon,'color','black')

figure
hold on
axesm miller
worldmap(latLim, lonLim)
geoshow(dlatx,dlonx,anom2_mean,'displaytype','texturemap');
load coastlines
geoshow(coastlat, coastlon)
colorbar
caxis([-30 30])
title('Total period Mean Temperature Amon')

 
 
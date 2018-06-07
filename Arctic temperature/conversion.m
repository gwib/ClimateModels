latLim = [58, 84];
lonLim = [-91, 29];

ncname1='./data/absolute_remap.nc'; % insert correct name
ncname2='./data/HadCRUT.4.6.0.0.median_remap.nc'; % insert correct name
ncdisp(ncname1);
finfo=ncinfo(ncname1);
dimNames = {finfo.Dimensions.Name};
varNames1 = {finfo.Variables.Name};
disp(dimNames);
disp(varNames1);
absTemp = ncread(ncname1,'tem');
lat = ncread(ncname1,'lat');
lon = ncread(ncname1,'lon');
dtime_hadCrut = ncread(ncname2, 'time');
%To plot a 2D matrix with geoshow the latitudes and longitudes must be converted to matrixes
for i=1:length(lat)
  dlonx_hadCrut(i,:)=lon;
end
dlonx_hadCrut = dlonx_hadCrut';
for i=1:length(lon)
  dlatx_hadCrut(:,i)=lat;
end
dlatx_hadCrut = dlatx_hadCrut';
ncdisp(ncname2);
finfo=ncinfo(ncname2);
dimNames = {finfo.Dimensions.Name};
varNames = {finfo.Variables.Name};
disp(dimNames);
disp(varNames);
anom = ncread(ncname2,'temperature_anomaly');
% m = repmat([1:12],1,168); m = cat(2,m,[1 2 3]);
abs2 = repmat(absTemp,1,1,168); %168 full years
abs2 = cat(3,abs2,absTemp(:,:,1),absTemp(:,:,2),absTemp(:,:,3));
for i = 1:2019
    hadCrutAbs(:,:,i) = anom(:,:,i)+abs2(:,:,i);
end
hadCrutAbs = hadCrutAbs(:,:,1561:1872);
hadCrutAbs_mean = nanmean(hadCrutAbs,3);

%% remove cells containing nan
nan_filter = double(isnan(hadCrutAbs));
nan_filter(nan_filter == 1) = nan;
nan_filter(nan_filter == 0) = 1;
 
% figure
% hold on
% worldmap world
% geoshow(dlatx_hadCrut,dlonx_hadCrut,hadCrutAbs_mean,'displaytype','texturemap');
% colorbar
% caxis([-30 30])
% load coastlines
% geoshow(coastlat,coastlon,'color','black')
%%
% figure
% hold on
% axesm miller
% worldmap(latLim, lonLim)
% geoimg = geoshow(dlatx_hadCrut,dlonx_hadCrut,hadCrutAbs_mean,'displaytype','texturemap');
% geoimg.AlphaDataMapping = 'none';
% geoimg.FaceAlpha = 'texturemap';
% alpha(geoimg,double(~isnan(hadCrutAbs_mean)))
% load coastlines
% geoshow(coastlat, coastlon)
% colorbar
% caxis([-50 20])
% title('Total period Mean Temperature HadCRUT')

 
 
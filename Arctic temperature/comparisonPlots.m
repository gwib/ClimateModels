% Plot total Temperature Mean
figure

subplot(2,2,1)
hold on
axesm miller
worldmap(latLim, lonLim)
geoimghadCrut = geoshow(dlatx_hadCrut,dlonx_hadCrut,hadCrutAbs_mean,'displaytype','texturemap');
geoimghadCrut.AlphaDataMapping = 'none';
geoimghadCrut.FaceAlpha = 'texturemap';
alpha(geoimghadCrut,double(~isnan(hadCrutAbs_mean)))
load coastlines
geoshow(coastlat, coastlon)
colorbar
caxis([-50 30])
title('Total period Mean Temperature HadCRUT')

subplot(2,2,2)
hold on
axesm miller
worldmap(latLim, lonLim)
geoshow(dlatx,dlonx,NCEP_obs_temp,'displaytype','texturemap');
load coastlines
geoshow(coastlat, coastlon)
colorbar
caxis([-50 30])
title('Total period Mean Temperature NCEP Reanalysis')

subplot(2,2,3)
hold on
axesm miller
worldmap(latLim, lonLim)
geoshow(dlatx,dlonx,TEMP_tas_Aamon_mean,'displaytype','texturemap');
load coastlines
geoshow(coastlat, coastlon)
colorbar
caxis([-50 30])
title('Total period Mean Temperature Amon')

subplot(2,2,4)
hold on
axesm miller
worldmap(latLim, lonLim)
geoimg_ArcHist = geoshow(dlatx,dlonx,tas_tasArcHistorical_mean,'displaytype','texturemap');
geoimg_ArcHist.AlphaDataMapping = 'none';
geoimg_ArcHist.FaceAlpha = 'texturemap';
alpha(geoimg_ArcHist,double(~isnan(tas_tasArcHistorical_mean)))
load coastlines
geoshow(coastlat, coastlon)
colorbar
caxis([-50 30])
title('Total period Mean Temperature Arc Historical (RCM)')

print('-dtiff','-r300','TotalMeanTemp');



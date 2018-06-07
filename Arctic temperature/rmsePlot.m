function rmsePlot(v1,v2, latLim, lonLim, dlatx, dlonx,rmseName)
   
    clim = [0 20];
    
    figure

    subplot(1,2,1)
    hold on
    axesm miller
    worldmap(latLim, lonLim)
    geoimghadCrut = geoshow(dlatx,dlonx,v1,'displaytype','texturemap');
    geoimghadCrut.AlphaDataMapping = 'none';
    geoimghadCrut.FaceAlpha = 'texturemap';
    alpha(geoimghadCrut,double(~isnan(v1)))
    load coastlines
    geoshow(coastlat, coastlon)
    %colormap(LB)
    caxis(clim)
    colorbar
    title(strcat('Reanalysis and GCM', rmseName))

    subplot(1,2,2)
    hold on
    axesm miller
    worldmap(latLim, lonLim)
    geoimghadCrut = geoshow(dlatx,dlonx,v2,'displaytype','texturemap');
    geoimghadCrut.AlphaDataMapping = 'none';
    geoimghadCrut.FaceAlpha = 'texturemap';
    alpha(geoimghadCrut,double(~isnan(v2)))
    load coastlines
    geoshow(coastlat, coastlon)
    %colormap(LB)
    caxis(clim)
    colorbar
    title(strcat('Reanalysis and RCM', rmseName))

    print('-dtiff','-r300',strcat('./plots/rmse_',rmseName));
end

function biasPlot(v1,v2,v3, biasName, latLim, lonLim, dlatx, dlonx)
   
    LB=flipud(lbmap(256,'BrownBlue'));
    m = max(abs(min(min(v1))),abs(max(max(v1))));
    clim = [-m m];
    
    figure

    subplot(2,2,1.5)
    hold on
    axesm miller
    worldmap(latLim, lonLim)
    geoimghadCrut = geoshow(dlatx,dlonx,v1,'displaytype','texturemap');
    geoimghadCrut.AlphaDataMapping = 'none';
    geoimghadCrut.FaceAlpha = 'texturemap';
    alpha(geoimghadCrut,double(~isnan(v1)))
    load coastlines
    geoshow(coastlat, coastlon)
    colormap(LB)
    caxis(clim)
    colorbar
    title('Observations and Reanalysis')

    subplot(2,2,3)
    hold on
    axesm miller
    worldmap(latLim, lonLim)
    geoimghadCrut = geoshow(dlatx,dlonx,v2,'displaytype','texturemap');
    geoimghadCrut.AlphaDataMapping = 'none';
    geoimghadCrut.FaceAlpha = 'texturemap';
    alpha(geoimghadCrut,double(~isnan(v2)))
    load coastlines
    geoshow(coastlat, coastlon)
    colormap(LB)
    caxis(clim)
    colorbar
    title('Reanalysis and GCM')

    subplot(2,2,4)
    hold on
    axesm miller
    worldmap(latLim, lonLim)
    geoimghadCrut = geoshow(dlatx,dlonx,v3,'displaytype','texturemap');
    geoimghadCrut.AlphaDataMapping = 'none';
    geoimghadCrut.FaceAlpha = 'texturemap';
    alpha(geoimghadCrut,double(~isnan(v3)))
    load coastlines
    geoshow(coastlat, coastlon)
    colormap(LB)
    caxis(clim)
    colorbar
    title('Reanalyis and RCM')

    print('-dtiff','-r300',strcat('./plots/bias_',biasName));
end

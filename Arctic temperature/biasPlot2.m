function biasPlot2(v1,v2, biasName, latLim, lonLim, dlatx, dlonx)
   
    LB=flipud(lbmap(256,'BrownBlue'));
    m = max(abs(min(min(v1))),abs(max(max(v1))));
    clim = [-m m];
    
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
    colormap(LB)
    caxis(clim)
    colorbar
    title(strcat('Reanalysis and GCM in period ', biasName))

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
    colormap(LB)
    caxis(clim)
    colorbar
    title(strcat('Reanalysis and RCM in period ', biasName))

    print('-dtiff','-r300',strcat('./plots/bias_',biasName));
end

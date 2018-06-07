function rmsePlot2(v1,v2,v3,v4, latLim, lonLim, dlatx, dlonx, rmse2Name)
   
    LB=flipud(lbmap(256,'BrownBlue'));
    %m = max(abs(min(min(v1))),abs(max(max(v1))));
    %clim = [-m m];
    
    figure

    subplot(2,2,1)
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
    %caxis(clim)
    colorbar
    title('Reanalysis and Observations Jan1')

    subplot(2,2,2)
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
    %caxis(clim)
    colorbar
    title('Reanalysis and Observations Jan2')
    
    subplot(2,2,3)
    hold on
    axesm miller
    worldmap(latLim, lonLim)
    geoimghadCrut = geoshow(dlatx,dlonx,v3,'displaytype','texturemap');
    geoimghadCrut.AlphaDataMapping = 'none';
    geoimghadCrut.FaceAlpha = 'texturemap';
    alpha(geoimghadCrut,double(~isnan(v3)))
    load coastlines
    geoshow(coastlat, coastlon)
    %colormap(LB)
    %caxis(clim)
    colorbar
    title('Reanalysis and Observations Jul1')
    
    subplot(2,2,4)
    hold on
    axesm miller
    worldmap(latLim, lonLim)
    geoimghadCrut = geoshow(dlatx,dlonx,v4,'displaytype','texturemap');
    geoimghadCrut.AlphaDataMapping = 'none';
    geoimghadCrut.FaceAlpha = 'texturemap';
    alpha(geoimghadCrut,double(~isnan(v4)))
    load coastlines
    geoshow(coastlat, coastlon)
    %colormap(LB)
    %caxis(clim)
    colorbar
    title('Reanalysis and Observations Jul2')
    
    print('-dtiff','-r300',strcat('./plots/rmse_',rmse2Name));
end

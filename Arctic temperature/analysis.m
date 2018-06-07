%%
run('loadData')
run('conversion')
run('comparisonPlots')
% time vector days since 1-11-1979
%% observation compared to reanalysis and models

% bias of total mean

% obs and Reanalysis
obsReanalysis_mean_bias = NCEP_mean - hadCrutAbs_mean;
% obs and GCM
obsGCM_mean_bias = TEMP_tas_Aamon_mean - hadCrutAbs_mean;
% obs and RCM
obsRCM_mean_bias = tas_tasArcHistorical_mean - hadCrutAbs_mean;

ncepGCM_mean_bias = TEMP_tas_Aamon_mean - NCEP_mean;
ncepRCM_mean_bias = tas_tasArcHistorical_mean - NCEP_mean;

biasPlot(obsReanalysis_mean_bias,ncepGCM_mean_bias,ncepRCM_mean_bias,'NCEPtotalMean', ...
    latLim, lonLim, dlatx, dlonx)

% comparison of RCM and GCM
GCM_RCM_mean_bias = TEMP_tas_Aamon_mean - tas_tasArcHistorical_mean;
m = max(abs(min(min(GCM_RCM_mean_bias))),abs(max(max(GCM_RCM_mean_bias))));
clim = [-m m];
LB=flipud(lbmap(256,'BrownBlue'));
figure
hold on
axesm miller
worldmap(latLim, lonLim)
geoimghadCrut = geoshow(dlatx,dlonx,GCM_RCM_mean_bias,'displaytype','texturemap');
geoimghadCrut.AlphaDataMapping = 'none';
geoimghadCrut.FaceAlpha = 'texturemap';
alpha(geoimghadCrut,double(~isnan(GCM_RCM_mean_bias)))
load coastlines
geoshow(coastlat, coastlon)
colormap(LB)
caxis(clim)
colorbar
title('GCM - RCM')
print('-dtiff','-r300','./plots/GCMminusRCM');

%% NCEP vs obs
m = max(abs(min(min(obsReanalysis_mean_bias))),abs(max(max(obsReanalysis_mean_bias))));
clim = [-m m];
LB=flipud(lbmap(256,'BrownBlue'));
figure
hold on
axesm miller
worldmap(latLim, lonLim)
geoimghadCrut = geoshow(dlatx,dlonx,obsReanalysis_mean_bias,'displaytype','texturemap');
geoimghadCrut.AlphaDataMapping = 'none';
geoimghadCrut.FaceAlpha = 'texturemap';
alpha(geoimghadCrut,double(~isnan(obsReanalysis_mean_bias)))
load coastlines
geoshow(coastlat, coastlon)
colormap(LB)
caxis(clim)
colorbar
title('NCEP vs HadCRUT (observations)')
print('-dtiff','-r300','./plots/NCEPvsObs');
%% examination of January values

% extract values for period 1 - 1980-1990
JANUARY_period1 = 1;
for i = 1:10
    temp_obs_Jan1(:,:,i) = hadCrutAbs(:,:,JANUARY_period1);
    temp_RCM_Jan1(:,:,i) = tas_tasArcHistorical(:,:,JANUARY_period1);
    temp_GCM_Jan1(:,:,i) = TEMP_tas_Aamon(:,:,JANUARY_period1);
    temp_NCEP_Jan1(:,:,i) = TEMP_obs_NCEP(:,:,JANUARY_period1);
    JANUARY_Jan1 = JANUARY_period1+12;
end

JANUARY_period2 = ((1995-1980)*12);
for i = 1:10
    temp_obs_Jan2(:,:,i) = hadCrutAbs(:,:,JANUARY_period2);
    temp_RCM_Jan2(:,:,i) = tas_tasArcHistorical(:,:,JANUARY_period2);
    temp_GCM_Jan2(:,:,i) = TEMP_tas_Aamon(:,:,JANUARY_period2);
    temp_NCEP_Jan2(:,:,i) = TEMP_obs_NCEP(:,:,JANUARY_period2);
    JANUARY_Jan2 = JANUARY_period2+12;
end

% extract values for period 1 - 1980-1990
JULY_period1 = 1;
for i = 1:10
    temp_obs_Jul1(:,:,i) = hadCrutAbs(:,:,JULY_period1);
    temp_RCM_Jul1(:,:,i) = tas_tasArcHistorical(:,:,JULY_period1);
    temp_GCM_Jul1(:,:,i) = TEMP_tas_Aamon(:,:,JULY_period1);
    temp_NCEP_Jul1(:,:,i) = TEMP_obs_NCEP(:,:,JULY_period1);
    JULY_period1 = JULY_period1+12;
end

% extract values for period 2 - 1995-2005
JULY_period2 = ((1995-1980)*12);
for i = 1:10
    temp_obs_Jul2(:,:,i) = hadCrutAbs(:,:,JULY_period2);
    temp_RCM_Jul2(:,:,i) = tas_tasArcHistorical(:,:,JULY_period2);
    temp_GCM_Jul2(:,:,i) = TEMP_tas_Aamon(:,:,JANUARY_period2);
    temp_NCEP_Jul2(:,:,i) = TEMP_obs_NCEP(:,:,JULY_period2);
    JULY_period2 = JULY_period2+12;
end

RMSE_obsRCM_Jan1 = NaN(320,160);
RMSE_obsGCM_Jan1 = NaN(320,160);
RMSE_obsNCEP_Jan1 = NaN(320,160);
RMSE_obsRCM_Jan2 = NaN(320,160);
RMSE_obsGCM_Jan2 = NaN(320,160);
RMSE_obsNCEP_Jan2 = NaN(320,160);
RMSE_obsRCM_Jul1 = NaN(320,160);
RMSE_obsGCM_Jul1 = NaN(320,160);
RMSE_obsNCEP_Jul1 = NaN(320,160);
RMSE_obsRCM_Jul2 = NaN(320,160);
RMSE_obsGCM_Jul2 = NaN(320,160);
RMSE_obsNCEP_Jul2 = NaN(320,160);


for i = 1:320
    for j=1:160
        % observations
        y_obsJan1 = temp_obs_Jan1(i,j,:);
        y_obsJan2 = temp_obs_Jan2(i,j,:);
        y_obsJul1 = temp_obs_Jul1(i,j,:);
        y_obsJul2 = temp_obs_Jul2(i,j,:);
        % RCM
        y_RCMJan1 = temp_RCM_Jan1(i,j,:);
        y_RCMJan2 = temp_RCM_Jan2(i,j,:);
        y_RCMJul1 = temp_RCM_Jul1(i,j,:);
        y_RCMJul2 = temp_RCM_Jul2(i,j,:);
        %GCM
        y_GCMJan1 = temp_GCM_Jan1(i,j,:);
        y_GCMJan2 = temp_GCM_Jan2(i,j,:);
        y_GCMJul1 = temp_GCM_Jul1(i,j,:);
        y_GCMJul2 = temp_GCM_Jul2(i,j,:);
        %NCEP
        y_NCEPJan1 = temp_NCEP_Jan1(i,j,:);
        y_NCEPJan2 = temp_NCEP_Jan2(i,j,:);
        y_NCEPJul1 = temp_NCEP_Jul1(i,j,:);
        y_NCEPJul2 = temp_NCEP_Jul2(i,j,:);
        
        
        RMSE_NCEPGCM_Jan1(i,j) = sqrt(nanmean((y_GCMJan1 - y_NCEPJan1).^2));
        RMSE_NCEPGCM_Jan2(i,j) = sqrt(nanmean((y_GCMJan2 - y_NCEPJan2).^2));
        RMSE_NCEPGCM_Jul1(i,j) = sqrt(nanmean((y_GCMJul1 - y_NCEPJul1).^2));
        RMSE_NCEPGCM_Jul2(i,j) = sqrt(nanmean((y_GCMJul2 - y_NCEPJul2).^2));
        
        RMSE_NCEPRCM_Jan1(i,j) = sqrt(nanmean((y_RCMJan1 - y_NCEPJan1).^2));
        RMSE_NCEPRCM_Jan2(i,j) = sqrt(nanmean((y_RCMJan2 - y_NCEPJan2).^2));
        RMSE_NCEPRCM_Jul1(i,j) = sqrt(nanmean((y_RCMJul1 - y_NCEPJul1).^2));
        RMSE_NCEPRCM_Jul2(i,j) = sqrt(nanmean((y_RCMJul2 - y_NCEPJul2).^2));
        
        RMSE_obsNCEP_Jan1(i,j) = sqrt(nanmean((y_obsJan1 - y_NCEPJan1).^2));
        RMSE_obsNCEP_Jan2(i,j) = sqrt(nanmean((y_obsJan2 - y_NCEPJan2).^2));
        RMSE_obsNCEP_Jul1(i,j) = sqrt(nanmean((y_obsJul1 - y_NCEPJul1).^2));
        RMSE_obsNCEP_Jul2(i,j) = sqrt(nanmean((y_NCEPJul2 - y_obsJul2).^2));
    end
end
rmsePlot2(RMSE_obsNCEP_Jan1, RMSE_obsNCEP_Jan2, RMSE_obsNCEP_Jul1, RMSE_obsNCEP_Jul2, latLim, lonLim, dlatx, dlonx,'obsReanalysis')
%%
rmsePlot(RMSE_NCEPGCM_Jan1,RMSE_NCEPRCM_Jan1, latLim, lonLim, dlatx, dlonx,'Jan1')
rmsePlot(RMSE_NCEPGCM_Jan2,RMSE_NCEPRCM_Jan2, latLim, lonLim, dlatx, dlonx,'Jan2')
rmsePlot(RMSE_NCEPGCM_Jul1,RMSE_NCEPRCM_Jul1, latLim, lonLim, dlatx, dlonx,'Jul1')
rmsePlot(RMSE_NCEPGCM_Jul2,RMSE_NCEPRCM_Jul2, latLim, lonLim, dlatx, dlonx,'Jul2')
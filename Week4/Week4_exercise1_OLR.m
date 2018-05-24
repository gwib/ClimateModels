% Extract CERES OLR

% Directory (correct to your working folder)
ncname='./../data/CERES_OLR_200003-201012.nc'; % insert correct name

%Display the metadata for the file (go though to familiarize with the format)
ncdisp(ncname);

%Get file information (and show dimension names and variable names)
finfo=ncinfo(ncname);
dimNames = {finfo.Dimensions.Name};
varNames = {finfo.Variables.Name};
disp(dimNames);
disp(varNames);

%Load the variables
dlon_OLR_CERES=ncread(ncname,'lon');
dlat_OLR_CERES=ncread(ncname,'lat');
dtime_OLR_CERES=ncread(ncname,'time');
OLR_CERES=ncread(ncname,'toa_lw_all_mon');

%%

% Extract ERAI OLR and Air Temp /T2M


ncname1='./../data/ERAI_OLR.nc'; % insert correct name
ncname2='./../data/ERAI_AIRTEMP.nc'; % insert correct name

%Display the metadata for the file (go though to familiarize with the format)
ncdisp(ncname1);
%Get file information (and show dimension names and variable names)
finfo=ncinfo(ncname1);
dimNames = {finfo.Dimensions.Name};
varNames = {finfo.Variables.Name};
disp(dimNames);
disp(varNames);
%Load the variables
dlon = ncread(ncname1,'longitude');
dlat = ncread(ncname1,'latitude');
dtime = ncread(ncname1,'time');
OLR_ERAI = ncread(ncname1,'ttr');


%Display the metadata for the file (go though to familiarize with the format)
ncdisp(ncname2);
%Get file information (and show dimension names and variable names)
finfo=ncinfo(ncname2);
dimNames = {finfo.Dimensions.Name};
varNames = {finfo.Variables.Name};
disp(dimNames);
disp(varNames);
%Load the variables
dlon = ncread(ncname2,'longitude');
dlat = ncread(ncname2,'latitude');
dtime = ncread(ncname2,'time');
T2M_ERAI = ncread(ncname2,'t2m');


%%

W = cos(pi/180*dlat)'; %weights (each cell 0-1)
WSUM = sum(W); % sum of weights
W = W/WSUM; % weights to equal 1
sum(W)

OLR_ERAI_lon_mean = nanmean(OLR_ERAI(:,:,:),1); % lon mean
OLR_ERAI_lon_mean = squeeze(OLR_ERAI_lon_mean(1,:,:)); % squeeze dimension
T2M_ERAI_lon_mean = nanmean(T2M_ERAI(:,:,:),1); % lon mean
T2M_ERAI_lon_mean = squeeze(T2M_ERAI_lon_mean(1,:,:)); % squeeze dimension

W_time = repmat(W,420,1)'; % repeat weigths over time

OLR_ERAI_lon_mean_W = W_time.*OLR_ERAI_lon_mean; % multiply weights with lon means
OLR_ERAI_lon_mean_W_SUM = sum(OLR_ERAI_lon_mean_W); % Summarize over latitudes
OLR_ERAI_lon_mean_W_SUM_NEG = OLR_ERAI_lon_mean_W_SUM*-1; % make negative for regression
T2M_ERAI_lon_mean_W = W_time.*T2M_ERAI_lon_mean; % multiply weights with lon means
T2M_ERAI_lon_mean_W_SUM = sum(T2M_ERAI_lon_mean_W); % Summarize over latitudes


time = [1:1:420]; % create simple time array for plot

% Find regression values
p = polyfit(T2M_ERAI_lon_mean_W_SUM,OLR_ERAI_lon_mean_W_SUM_NEG,1);
REG_line_x = [min(T2M_ERAI_lon_mean_W_SUM) max(T2M_ERAI_lon_mean_W_SUM)]; % Data points for regresion line
REG_line_y = REG_line_x*p(1)+p(2);
Gamma = p(1);


figure
hold on
plot(time,OLR_ERAI_lon_mean_W_SUM_NEG);
ylabel('OLR (w/m2)');
xlabel('Time (months from beginning)');
print('-dtiff','-r300','OLR_plot');

figure
hold on
plot(time,T2M_ERAI_lon_mean_W_SUM);
ylabel('T2M (K)');
xlabel('Time (months from beginning)');
print('-dtiff','-r300','T2M_plot');

figure
hold on
plot(REG_line_x,REG_line_y,'r','linewidth',3);
plot(T2M_ERAI_lon_mean_W_SUM,OLR_ERAI_lon_mean_W_SUM_NEG,'*k');
xlabel('T2m');
ylabel('Negative OLR');
str=sprintf('Gamma = %.4f', p(1));
title(str)
print('-dtiff','-r300','Gamma_plot');

%%
% Same with ERAI T2M and CERES OLR

% Input for regression
% T2M_ERAI_lon_mean_W_SUM already given from script


OLR_CERES_lon_mean = nanmean(OLR_CERES(:,:,:),1); % lon mean
OLR_CERES_lon_mean = squeeze(OLR_CERES_lon_mean(1,:,:)); % squeeze dimension

W_CERES = cos(pi/180*dlat_OLR_CERES)'; %weights (each cell 0-1)
WSUM_CERES = sum(W_CERES); % sum of weights
W_CERES = W_CERES/WSUM_CERES; % weights to equal 1
W_time_CERES = repmat(W_CERES,130,1)';

OLR_CERES_lon_mean_W = W_time_CERES.*OLR_CERES_lon_mean; % multiply weights with lon means
OLR_CERES_lon_mean_W_SUM = sum(OLR_CERES_lon_mean_W); % Summarize over latitudes
OLR_CERES_lon_mean_W_SUM_NEG = OLR_CERES_lon_mean_W_SUM*-1; % make negative for regression

% Find regression values
%TODO: find overlapping time period
% CERES time period: 03/2000 - 12/2010
% ERAI time period: 12/1979 - 12/2014
% index of overlapping time period: 267:396


p2 = polyfit(T2M_ERAI_lon_mean_W_SUM(267:396),OLR_CERES_lon_mean_W_SUM_NEG,1);
REG_line_x2 = [min(T2M_ERAI_lon_mean_W_SUM(267:396)) max(T2M_ERAI_lon_mean_W_SUM(267:396))]; % Data points for regresion line
REG_line_y2 = REG_line_x2*p2(1)+p2(2);
Gamma = p2(1);

figure
hold on
plot(dtime_OLR_CERES,OLR_CERES_lon_mean_W_SUM_NEG);
ylabel('OLR (w/m2)');
xlabel('Time (months from beginning)');
print('-dtiff','-r300','OLR_plot');
title('OLR CERES')

figure
hold on
plot(dtime(267:396),T2M_ERAI_lon_mean_W_SUM(267:396));
ylabel('T2M (K)');
xlabel('Time (months from beginning)');
title('OLR CERES')
print('-dtiff','-r300','T2M_plot');

figure
hold on
erai_reg = plot(REG_line_x,REG_line_y,'r','linewidth',3);
plot(T2M_ERAI_lon_mean_W_SUM(267:396),OLR_ERAI_lon_mean_W_SUM_NEG(267:396),'*k');
erai_ceres_reg = plot(REG_line_x2,REG_line_y2,'g','linewidth',3);
plot(T2M_ERAI_lon_mean_W_SUM(267:396),OLR_CERES_lon_mean_W_SUM_NEG,'*k');

str1=sprintf('Gamma = %.4f', p(1));
label(erai_reg,str1)
str2=sprintf('Gamma = %.4f', p2(1));
label(erai_ceres_reg, str2)

xlabel('T2m');
ylabel('Negative OLR');

diffstr = fprintf('Difference in Gamma: %.4f',gammaDiff);
title(diffstr)
print('-dtiff','-r300','Gamma_plot');

gammaDiff = p2(1)-p(1);



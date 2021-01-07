%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% ADS-B TRUST - PLOTS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;
format long;

E = referenceEllipsoid('Earth');

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% IMPORT RESULTS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf("Importing Results ...\n");
tic

% Filepath
day = "2020-02-15";
%filedir_reports = fullfile("../Reports", day);
filedir_Data = fullfile("../Data");
%filedir_Grid = fullfile("../Grid");
filedir_Results = fullfile("../Results");
filedir_Figures = fullfile("Figures");

if ~exist(filedir_Figures, "dir")
    fprintf("\nGenerating Directory for Figures ...\n");
    
    % Create Figure Folder
    mkdir("Figures");
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% TOTAL MESSAGES %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
load(fullfile(filedir_Data, "Grid.mat"))

% Figure
figure
title(" ")
worldmap("Europe")

% Calculate Borders
for lat = 1:numel(grid_lat)
    tmp_lat = []; tmp_lon = [];
    
    % Latitude Borders
    [tmp_lat(1),~] = reckon('rh', grid_lat(lat), 0, 0.5 * grid_resolution, 0, E);
    [tmp_lat(2),~] = reckon('rh', grid_lat(lat), 0, 0.5 * grid_resolution, 180, E);
    
    % Longitude Borders
    [~,tmp_lon(1)] = reckon('rh', grid_lat(lat), grid_lon{lat}(1), 0.5 * grid_resolution, -90, E);
    
    for lon = 2:numel(grid_lon{lat})
        [~,tmp(1)] = reckon('rh', grid_lat(lat), grid_lon{lat}(lon-1), 0.5 * grid_resolution, 90, E);
        [~,tmp(2)] = reckon('rh', grid_lat(lat), grid_lon{lat}(lon), 0.5 * grid_resolution, -90, E);
    
        tmp_lon(lon) = 0.5 * (tmp(1)+tmp(2));
    end
    [~,tmp_lon(numel(grid_lon{lat})+1)] = reckon('rh', grid_lat(lat), grid_lon{lat}(end), 0.5 * grid_resolution, 90, E);
    
    surfacem(tmp_lat,tmp_lon,log10(double(Grid(lat).Total_Messages)+1), 'LineStyle','none');
end

% Coastlines
load coastlines
plotm(coastlat,coastlon, 'Color','black', 'LineWidth',1);

% Color
j = jet;
j(1,:) = [1 1 1];
colormap(j);
cbar = colorbar('northoutside');
caxis([0 5])
cbar.TickLabels = {'0','10^1','10^2','10^3','10^4','>10^5'};
ylabel(cbar, "Number of Messages")

% Save Figure
print('-opengl','-dpdf',fullfile(filedir_Figures, "Total_Messages_color"),'-r300')

% Gray
j = gray;
j = j(end:-1:1,:);
colormap(j);
cbar = colorbar('northoutside');
caxis([0 5])
cbar.TickLabels = {'0','10^1','10^2','10^3','10^4','>10^5'};
ylabel(cbar, "Number of Messages")
print('-opengl','-dpdf',fullfile(filedir_Figures, "Total_Messages_gray"),'-r300')

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% SENSOR COVERAGE %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
title(' ');
worldmap europe

% Calculate Borders
for lat = 1:numel(grid_lat)
    tmp_lat = []; tmp_lon = [];
    
    % Latitude Borders
    [tmp_lat(1),~] = reckon('rh', grid_lat(lat), 0, 0.5 * grid_resolution, 0, E);
    [tmp_lat(2),~] = reckon('rh', grid_lat(lat), 0, 0.5 * grid_resolution, 180, E);
    
    % Longitude Borders
    [~,tmp_lon(1)] = reckon('rh', grid_lat(lat), grid_lon{lat}(1), 0.5 * grid_resolution, -90, E);
    
    for lon = 2:numel(grid_lon{lat})
        [~,tmp(1)] = reckon('rh', grid_lat(lat), grid_lon{lat}(lon-1), 0.5 * grid_resolution, 90, E);
        [~,tmp(2)] = reckon('rh', grid_lat(lat), grid_lon{lat}(lon), 0.5 * grid_resolution, -90, E);
    
        tmp_lon(lon) = 0.5 * (tmp(1)+tmp(2));
    end
    [~,tmp_lon(numel(grid_lon{lat})+1)] = reckon('rh', grid_lat(lat), grid_lon{lat}(end), 0.5 * grid_resolution, 90, E);
    
    surfacem(tmp_lat,tmp_lon,sum(logical(Grid(lat).Sensors)), 'LineStyle','none');
end

% Coastlines
load coastlines
plotm(coastlat,coastlon, 'Color','black', 'LineWidth',1);

% Color
j = jet;
j(1,:) = [1 1 1];
colormap(j);
cbar = colorbar('northoutside');
caxis([0 10])
cbar.TickLabels = {'0','2','4','6','8','>10'};
ylabel(cbar, "Sensor Coverage")

% Save Figure
print('-opengl','-dpdf',fullfile(filedir_Figures, "Sensor_Coverage_color"),'-r300')

% Gray
j = gray;
j = j(end:-1:1,:);
colormap(j);
cbar = colorbar('northoutside');
caxis([0 10])
cbar.TickLabels = {'0','2','4','6','8','>10'};
ylabel(cbar, "Sensor Coverage")
print('-opengl','-dpdf',fullfile(filedir_Figures, "Sensor_Coverage_gray"),'-r300')

%%%%%%%%%%%%%%%%%%%%%%
%%%%% RESILIENCE %%%%%
%%%%%%%%%%%%%%%%%%%%%%
figure
title(' ');
worldmap europe

% Calculate Borders
for lat = 1:numel(grid_lat)
    tmp_lat = []; tmp_lon = [];
    
    % Latitude Borders
    [tmp_lat(1),~] = reckon('rh', grid_lat(lat), 0, 0.5 * grid_resolution, 0, E);
    [tmp_lat(2),~] = reckon('rh', grid_lat(lat), 0, 0.5 * grid_resolution, 180, E);
    
    % Longitude Borders
    [~,tmp_lon(1)] = reckon('rh', grid_lat(lat), grid_lon{lat}(1), 0.5 * grid_resolution, -90, E);
    
    for lon = 2:numel(grid_lon{lat})
        [~,tmp(1)] = reckon('rh', grid_lat(lat), grid_lon{lat}(lon-1), 0.5 * grid_resolution, 90, E);
        [~,tmp(2)] = reckon('rh', grid_lat(lat), grid_lon{lat}(lon), 0.5 * grid_resolution, -90, E);
    
        tmp_lon(lon) = 0.5 * (tmp(1)+tmp(2));
    end
    [~,tmp_lon(numel(grid_lon{lat})+1)] = reckon('rh', grid_lat(lat), grid_lon{lat}(end), 0.5 * grid_resolution, 90, E);
    
    surfacem(tmp_lat,tmp_lon,sum(logical(Grid(lat).Sensors)), 'LineStyle','none');
end

% Coastlines
load coastlines
plotm(coastlat,coastlon, 'Color','black', 'LineWidth',1);

% Color
j = jet;
colormap(flipud(j(129:end,:)));
cbar = colorbar('northoutside');
caxis([0 10])
cbar.TickLabels = {'0','2','4','6','8','>10'};

% Save Figure
print('-opengl','-dpdf',fullfile(filedir_Figures, "Resilience_color"),'-r300')
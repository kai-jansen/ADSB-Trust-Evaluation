%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% ADS-B TRUST - 10.11.2020 %%%%%
%%%%%%%%%%%% Kai Jansen %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;
format long;

E = referenceEllipsoid('Earth');

%%%%%%%%%%%%%%%%%%%%%%%
%%%%% IMPORT DATA %%%%%
%%%%%%%%%%%%%%%%%%%%%%%
fprintf("Importing Reports ...\n");
tic

% Filepath
%day = "2020-02-15";
day = "2020-02-17";

% Internal Storage
%filedir_Reports = fullfile("..", day, "Reports");
filedir_Data = fullfile("..", day, "Data");
%filedir_Grid = fullfile("..", day, "Grid");
filedir_Results = fullfile("..", day, "Results");

% External Storage


% Data Format of *.csv Files
report_format = {'%u','%6c','%f','%f','%f','%f','%f','%s','%C','%C','%C','%u16','%f','%f','%f','%f','%s'};

% Datastore
ds = tabularTextDatastore(fullfile(filedir_Reports, "*.csv"));
ds.SelectedFormats = report_format;
ds.SelectedVariableNames = {'time','icao24','lat','lon','sensors'}; %,'geoaltitude','baroaltitude'
T = tall(ds);

% Number of Reports
%heightT = gather(height(T));
%heightT = 132883464; % 2020-02-15
heightT = 134918966; % 2020-02-17

fprintf("Reference to Reports successfully set! [Time: %s]\n", seconds(toc));

%%%%%%%%%%%%%%%%%%%%%&%%%%%
%%%%% SENSOR ANALYSIS %%%%%
%%%%%%%%%%%%%%%%%%%%%%&%%%%
if ~exist(fullfile(filedir_Data, "Sensors.mat"), "file")
    fprintf("\nDetermining the Number of Sensors ...\n");
    
    % Sensors
    Sensors = [];
    process_at_once = 1000000;
    
    for hour = 1:24
        tic
        
        % Datastore
        ds = tabularTextDatastore(fullfile(filedir_Reports, num2str(hour-1,'%02d') + ".csv"));
        ds.SelectedFormats = report_format;
        ds.SelectedVariableNames = {'sensors'};
        D = tall(ds);
        
        % Read Receiving Sensors
        Sensors_received = gather(D.sensors);
        
        % Process Data in Chunks
        reports_processed = 0;
        while reports_processed < numel(Sensors_received)
            % Find Unique Sensors
            if reports_processed + process_at_once <= numel(Sensors_received)
                Sensors = unique([Sensors;split(join(Sensors_received(reports_processed+1:reports_processed+process_at_once),';'),';')]);
                reports_processed = reports_processed + process_at_once;
            else
                Sensors = unique([Sensors;split(join(Sensors_received(reports_processed+1:end),';'),';')]);
                reports_processed = numel(Sensors_received);
            end
        end
        fprintf("Processed: %d out of 24 Files [Time: %s | Remaining: %s]\n\n", ...
                hour, seconds(toc), duration(0,0,toc*(24-hour)));
    end
    save(fullfile(filedir_Data, "Sensors.mat"), "Sensors")
else
    load(fullfile(filedir_Data, "Sensors.mat"))
end
fprintf("Sensors successfully analyzed! [Time: %s]\n", seconds(toc));

%%%%%%%%%%%%%%%%%%%%%%%
%%%%% EUROPE GRID %%%%%
%%%%%%%%%%%%%%%%%%%%%%%
if ~exist(fullfile(filedir_Data, "Grid.mat"), "file")
    fprintf("\nGenerating Grid ...\n");
    tic
    
    % Borders
    lat_min = 30; lat_max = 75;
    lon_min = -25; lon_max = 45;
    grid_resolution = 10e3; % 10 [km]
    
    % Grid Lat
    grid_lat_dist = distance('rh', lat_max, 0, lat_min, 0, E);
    grid_lat = NaN(ceil(grid_lat_dist / grid_resolution),1);
    
    % Equally Divide Lat
    lat_out_of_border = 0.5 * (numel(grid_lat) * grid_resolution - grid_lat_dist);
    [grid_lat(1),~] = reckon('rh', lat_max, 0, 0.5 * grid_resolution - lat_out_of_border, 180, E);
    for i = 2:numel(grid_lat)
        [grid_lat(i),~] = reckon('rh', grid_lat(i-1), 0, grid_resolution, 180, E);
    end
    
    % Grid Lon
    grid_lon = cell(numel(grid_lat),1);
    
    % For Each Lat Index
    for i = 1:numel(grid_lat)
        grid_lon_dist = distance('rh', grid_lat(i), lon_min, grid_lat(i), lon_max, E);
        grid_lon_tmp = NaN(1,ceil(grid_lon_dist / grid_resolution));
        
        % Equally Divide Lon
        lon_out_of_border = 0.5 * (numel(grid_lon_tmp) * grid_resolution - grid_lon_dist);
        [~,grid_lon_tmp(1)] = reckon('rh', grid_lat(i), lon_min, 0.5 * grid_resolution - lon_out_of_border, 90, E);
        
        for j = 2:numel(grid_lon_tmp)
            [~,grid_lon_tmp(j)] = reckon('rh', grid_lat(i), grid_lon_tmp(j-1), grid_resolution, 90, E);
        end
        
        % Save in Cell Array
        grid_lon{i} = grid_lon_tmp;
    end
    
    % Grid
    for i = 1:numel(grid_lat)
        Grid(i).lat = grid_lat(i);
        Grid(i).lon = grid_lon{i};
        Grid(i).Sensors = zeros(numel(Sensors),numel(grid_lon{i}),'uint32');
        Grid(i).Total_Messages = zeros(1,numel(grid_lon{i}),'uint32');
    end
    
    save(fullfile(filedir_Data, "Grid.mat"), "Grid", "grid_lat", "grid_lon", "grid_resolution")
else
    load(fullfile(filedir_Data, "Grid.mat"))
end
fprintf("Grid successfully generated! [Time: %s]\n", seconds(toc));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INDEX COMPUTATION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\nComputing Indices ...\n');

% Find all unique Aircraft
if ~exist(fullfile(filedir_Data, "index_icao.mat"), "file")
    tic
    
    index_icao = zeros(heightT,1,'uint16');
    tmp_icao = gather(T.icao24);
    aircraft_unique = unique(tmp_icao, 'rows');
    
    % ICAO Index
    [~,index_icao(:)] = ismember(tmp_icao, aircraft_unique, 'rows');
    clear tmp_icao
    
    save(fullfile(filedir_Data, "index_icao.mat"), "index_icao", "aircraft_unique", '-v7.3')
    fprintf("Index for aircraft successfully computed! [Time: %s]\n\n", seconds(toc));
else
    load(fullfile(filedir_Data, "index_icao.mat"))
end

% Calculate Lat Index
if ~exist(fullfile(filedir_Data, "index_lat.mat"), "file")
    tic
    
    index_lat = zeros(heightT,1,'uint16');
    tmp_lat = gather(T.lat);
    
    parfor i = 1:heightT
        % Distance to all Lat in Grid
        [~,index_lat(i)] = min(distance('rh', tmp_lat(i), 0, grid_lat, 0, E));
    end
    clear tmp_lat
    
    save(fullfile(filedir_Data, "index_lat.mat"), "index_lat", '-v7.3')
    fprintf("Index for latitude successfully computed! [Time: %s]\n\n", seconds(toc));
else
    load(fullfile(filedir_Data, "index_lat.mat"))
end

% Calculate Lon Index
if ~exist(fullfile(filedir_Data, "index_lon.mat"), "file")
    tic
    
    index_lon = zeros(heightT,1,'uint16');
    tmp_lon = gather(T.lon);
    
    parfor i = 1:heightT
        % Distance to all Lon in Grid
        [~,index_lon(i)] = min(distance('rh', Grid(index_lat(i)).lat, tmp_lon(i), ...
                                              Grid(index_lat(i)).lat, Grid(index_lat(i)).lon, E));
    end
    clear tmp_lon
    
    save(fullfile(filedir_Data, "index_lon.mat"), "index_lon", '-v7.3')
    fprintf("Index for longitude successfully computed! [Time: %s]\n\n", seconds(toc));
else
    load(fullfile(filedir_Data, "index_lon.mat"))
end

% Compute Indicies for Sensors
if ~exist(fullfile(filedir_Data, "index_sensors.mat"), "file")
    tic
    
    index_sensors = zeros(heightT,ceil(numel(Sensors)/8),'uint8');
    
    % LUT for Efficient Bit Setting
    LUT_ADD = zeros(numel(Sensors),size(index_sensors,2),'uint8');
    for i = 1:size(index_sensors,2)
        for j = 1:8
            LUT_ADD((i-1)*8+j,i) = bitset(LUT_ADD((i-1)*8+j,i),j);
            
            % Catch End
            if (i-1)*8+j == numel(Sensors)
                break
            end
        end
    end
    
    % Process all Reports
    reports_processed = 0;
    for hour = 1:24
        tic
        
        % Datastore
        ds = tabularTextDatastore(fullfile(filedir_Reports, num2str(hour-1,'%02d') + ".csv"));
        ds.SelectedFormats = report_format;
        ds.SelectedVariableNames = {'sensors'};
        D = tall(ds);
        
        % Read Receiving Sensors
        Sensors_received = gather(D.sensors);
        
        for i = 1:numel(Sensors_received)
            [~,idx] = ismember(split(Sensors_received(i),';'), Sensors);
            
            % Sensors Index
            index_sensors(reports_processed+i,:) = sum(LUT_ADD(idx,:),1,'native');
            
            % Grid Update
            Grid(index_lat(reports_processed+i)).Sensors(idx,index_lon(reports_processed+i)) = ...
                Grid(index_lat(reports_processed+i)).Sensors(idx,index_lon(reports_processed+i)) + 1;
            Grid(index_lat(reports_processed+i)).Total_Messages(index_lon(reports_processed+i)) = ...
                Grid(index_lat(reports_processed+i)).Total_Messages(index_lon(reports_processed+i)) + 1;
        end
        
        % Processed Reports
        reports_processed = reports_processed + numel(Sensors_received);
        
        fprintf("Processed: %d out of 24 Files [Time: %s | Remaining: %s]\n\n", ...
                hour, seconds(toc), duration(0,0,toc*(24-hour)));
    end
    save(fullfile(filedir_Data, "index_sensors.mat"), "index_sensors", '-v7.3')
    save(fullfile(filedir_Data, "Grid.mat"), "Grid", "grid_lat", "grid_lon", "grid_resolution")
    fprintf("Index for sensors successfully computed! [Time: %s]\n", seconds(toc));
else
    load(fullfile(filedir_Data, "index_sensors.mat"))
end
fprintf("Indices for aircraft, grid, and sensors successfully computed! [Time: %s]\n", seconds(toc));

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% REPORT MAPPING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist(filedir_Grid, "dir")
    fprintf("\nMapping Reports to Grid ...\n");
    
    % Create Grid Folder
    mkdir(filedir_Grid);
    
    % Datastore
    ds = tabularTextDatastore(fullfile(filedir_Reports, "*.csv"));
    ds.SelectedFormats = report_format;
    ds.SelectedVariableNames = {'time','lat','lon'};
    D = tall(ds);
    
    % Read Time Latitude Longitude
    [tmp_time,tmp_lat,tmp_lon] = gather(D.time, D.lat, D.lon);
    
    % Process all Reports
    reports_processed = 0;
    
    % For all Lat Indices
    for lat = 1:numel(grid_lat)
        tic
        
        % Reports in Grid Lat Index
        select_lat = index_lat==lat;
        
        % Check if any Reports in Grid Lat Index
        if ~any(select_lat)
            continue
        end
        
        % Make New Folder for Lat Index
        mkdir(filedir_Grid, num2str(lat));
        
        % For all Lon Indices
        for lon = 1:numel(grid_lon{lat})
            
            % Reports in Grid Lon Index
            select = (select_lat & index_lon==lon);

            % Check if any Reports in Grid Lon Index
            if ~any(select)
                continue
            end

            % All Reports in Grid Index
            Normal_Operation = table(tmp_time(select), ...
                                     index_icao(select), ...
                                     tmp_lat(select), ...
                                     tmp_lon(select), ...
                                     index_sensors(select,:), ...
                                     'VariableNames', {'time','icao24','lat','lon','sensors'});

            % Processed Reports
            reports_processed = reports_processed + sum(select);

            % Make New Folder for Lon Index
            mkdir(fullfile(filedir_Grid, num2str(lat)), num2str(lon));

            % Save Data
            writetable(Normal_Operation, fullfile(filedir_Grid, num2str(lat), num2str(lon), "Normal_Operation.csv"));
        end
        fprintf("Latitude Index %d out of %d (%.2f%%) successfully mapped! [Time: %s | Remaining: %s]\n", ...
                lat, numel(grid_lat), 100 * lat / numel(grid_lat), ...
                seconds(toc), duration(0,0,toc*(numel(grid_lat)-lat)));
    end
    fprintf("Reports successfully mapped to Grid!\n");
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% SIMULATION: GPS SPOOFING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf("\nGenerating ADS-B Reports for Attack: GPS Spoofing ...\n");

% Datastore
ds = tabularTextDatastore(fullfile(filedir_Reports, "*.csv"));
ds.SelectedFormats = report_format;
ds.SelectedVariableNames = {'time','lat','lon'};
D = tall(ds);

% Read Time Latitude Longitude
[tmp_time,tmp_lat,tmp_lon] = gather(D.time, D.lat, D.lon);

% Generate at least the Number of ADS-B Reports for Normal Operation
GPS_Spoofing_num = heightT;
GPS_Spoofing_done = 0;

% GPS Spoofing Data
GPS_Spoofing_idx = zeros(heightT,1,'uint32');
GPS_Spoofing_icao = zeros(heightT,1,'uint16');
GPS_Spoofing_lat = zeros(heightT,1);
GPS_Spoofing_lon = zeros(heightT,1);
GPS_Spoofing_index_lat = zeros(heightT,1,'uint16');
GPS_Spoofing_index_lon = zeros(heightT,1,'uint16');
GPS_Spoofing_deviation = zeros(heightT,1,'uint8');
GPS_Spoofing_duration = zeros(heightT,1,'uint16');
GPS_Spoofing_distance = zeros(heightT,1);

while GPS_Spoofing_done < GPS_Spoofing_num
    tic
    
    % Random Aircraft Report
    select = randi(heightT,1);
    
    % Indices for Reports within the Next Hour
    GPS_Spoofing_idx_tmp = find(index_icao == index_icao(select) & ...
                                tmp_time >= tmp_time(select) & ...
                                tmp_time <= (tmp_time(select) + 3600));
    
    % Reports for Track
    track = table(tmp_time(GPS_Spoofing_idx_tmp), ...
                  index_icao(GPS_Spoofing_idx_tmp), ...
                  tmp_lat(GPS_Spoofing_idx_tmp), ...
                  tmp_lon(GPS_Spoofing_idx_tmp), ...
                  index_sensors(GPS_Spoofing_idx_tmp,:), ...
                  'VariableNames', {'time','icao24','lat','lon','sensors'});
    
    % Duration
    GPS_Spoofing_duration_tmp = track.time - track.time(1);
    GPS_Spoofing_duration_tmp = repmat(GPS_Spoofing_duration_tmp,2,1);
    
    % Distances and Angles to all the next Points
    [arclen,az] = distance(track.lat(1), track.lon(1), ...
                           track.lat, track.lon, E);
    
    % For all Deviations
    for deviation = [1, 2, 5, 10, 20, 45]
        
        % Calculate new Positions with Deviation of x°
        [tmp_lat_plus,tmp_lon_plus] = reckon(track.lat(1), track.lon(1), arclen, az+deviation, E);
        [tmp_lat_minus,tmp_lon_minus] = reckon(track.lat(1), track.lon(1), arclen, az-deviation, E);
        GPS_Spoofing_lat_tmp = [tmp_lat_plus; tmp_lat_minus];
        GPS_Spoofing_lon_tmp = [tmp_lon_plus; tmp_lon_minus];
        
        % Distances
        GPS_Spoofing_distance_tmp = distance(track.lat, track.lon, ...
                                             tmp_lat_plus, tmp_lon_plus, E);
        GPS_Spoofing_distance_tmp = repmat(GPS_Spoofing_distance_tmp,2,1);
        
        % Calculate Grid Index
        GPS_Spoofing_index_lat_tmp = zeros(numel(GPS_Spoofing_lat_tmp),1,'uint16');
        GPS_Spoofing_index_lon_tmp = zeros(numel(GPS_Spoofing_lon_tmp),1,'uint16');
        for i = 1:numel(GPS_Spoofing_lat_tmp)
            % Distance to all Lat in Grid
            [~,GPS_Spoofing_index_lat_tmp(i)] = min(distance('rh', GPS_Spoofing_lat_tmp(i), 0, ...
                                                                   grid_lat, 0, E));
            % Distance to all Lon in Grid
            [~,GPS_Spoofing_index_lon_tmp(i)] = min(distance('rh', Grid(GPS_Spoofing_index_lat_tmp(i)).lat, GPS_Spoofing_lon_tmp(i), ...
                                                                   Grid(GPS_Spoofing_index_lat_tmp(i)).lat, Grid(GPS_Spoofing_index_lat_tmp(i)).lon, E));
        end
        
        % Check if Area is within Grid and Different
        idx = GPS_Spoofing_index_lat_tmp ~= 0 & ...
              GPS_Spoofing_index_lon_tmp ~= 0 & ...
             (GPS_Spoofing_index_lat_tmp ~= repmat(index_lat(GPS_Spoofing_idx_tmp),2,1) | ...
              GPS_Spoofing_index_lon_tmp ~= repmat(index_lon(GPS_Spoofing_idx_tmp),2,1));
        
        % Discard Spoofing Reports in Grids without authentic Reports
        for i = 1:numel(idx)
            if idx(i)
                if ~Grid(GPS_Spoofing_index_lat_tmp(i)).Total_Messages(GPS_Spoofing_index_lon_tmp(i))
                    idx(i) = false;
                end
            end
        end
        
        % GPS Spoofing Reports
        tmp_idx = repmat(GPS_Spoofing_idx_tmp,2,1);
        GPS_Spoofing_idx(GPS_Spoofing_done+1:GPS_Spoofing_done+nnz(idx)) = tmp_idx(idx);
        
        GPS_Spoofing_icao(GPS_Spoofing_done+1:GPS_Spoofing_done+nnz(idx)) = repmat(index_icao(select),nnz(idx),1);
        GPS_Spoofing_lat(GPS_Spoofing_done+1:GPS_Spoofing_done+nnz(idx)) = GPS_Spoofing_lat_tmp(idx);
        GPS_Spoofing_lon(GPS_Spoofing_done+1:GPS_Spoofing_done+nnz(idx)) = GPS_Spoofing_lon_tmp(idx);
        GPS_Spoofing_index_lat(GPS_Spoofing_done+1:GPS_Spoofing_done+nnz(idx)) = GPS_Spoofing_index_lat_tmp(idx);
        GPS_Spoofing_index_lon(GPS_Spoofing_done+1:GPS_Spoofing_done+nnz(idx)) = GPS_Spoofing_index_lon_tmp(idx);
        
        GPS_Spoofing_deviation(GPS_Spoofing_done+1:GPS_Spoofing_done+nnz(idx)) = repmat(deviation,nnz(idx),1);
        GPS_Spoofing_duration(GPS_Spoofing_done+1:GPS_Spoofing_done+nnz(idx)) = GPS_Spoofing_duration_tmp(idx);
        GPS_Spoofing_distance(GPS_Spoofing_done+1:GPS_Spoofing_done+nnz(idx)) = GPS_Spoofing_distance_tmp(idx);
        
        % Reports Generated
        GPS_Spoofing_done = GPS_Spoofing_done + nnz(idx);
    end 
    fprintf("Successfully generated GPS Spoofing reports! %d out of %d (%.2f%%) [Time: %s]\n", ...
            GPS_Spoofing_done, heightT, 100 * GPS_Spoofing_done / heightT, seconds(toc));
end
clear GPS_Spoofing_idx_tmp
clear GPS_Spoofing_lat_tmp
clear GPS_Spoofing_lon_tmp
clear GPS_Spoofing_index_lat_tmp
clear GPS_Spoofing_index_lon_tmp
clear GPS_Spoofing_duration_tmp
clear GPS_Spoofing_distance_tmp

% Unique Grid Indices
[grid_index_unique,~,idx_tmp] = unique([GPS_Spoofing_index_lat, GPS_Spoofing_index_lon], 'rows');

GPS_Spoofing_done = 0;
% Save for each Unique Index
for i = 1:size(grid_index_unique,1)
    tic
    
    tmp = find(idx_tmp==i);
    % Table for Spoofing Reports
    GPS_Spoofing = table(tmp_time(GPS_Spoofing_idx(tmp)), ...
                         GPS_Spoofing_icao(tmp), ...
                         GPS_Spoofing_lat(tmp), ...
                         GPS_Spoofing_lon(tmp), ...
                         index_sensors(GPS_Spoofing_idx(tmp),:), ...
                         GPS_Spoofing_deviation(tmp), ...
                         GPS_Spoofing_duration(tmp), ...
                         GPS_Spoofing_distance(tmp), ...
                         'VariableNames', {'time','icao24','lat','lon','sensors','deviation','duration','distance'});
    
    % Save Data
    writetable(GPS_Spoofing, fullfile(filedir_Grid, num2str(grid_index_unique(i,1)), num2str(grid_index_unique(i,2)), "GPS_Spoofing.csv"));
    
    % Reports Saved
    GPS_Spoofing_done = GPS_Spoofing_done + numel(tmp);
    
    fprintf("Successfully saved GPS Spoofing reports! %d out of %d (%.2f%%) [Time: %s | Remaining: %s]\n", ...
            GPS_Spoofing_done, size(GPS_Spoofing_idx,1), 100 * GPS_Spoofing_done / size(GPS_Spoofing_idx,1), ...
            seconds(toc), duration(0,0,toc*(size(grid_index_unique,1)-i)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% SIMULATION: ADS-B SPOOFING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf("\nGenerating ADS-B Reports for Attack: ADS-B Spoofing ...\n");

% Datastore
ds = tabularTextDatastore(fullfile(filedir_Reports, "*.csv"));
ds.SelectedFormats = report_format;
ds.SelectedVariableNames = {'time','lat','lon'};
D = tall(ds);

% Read Time Latitude Longitude
[tmp_time,tmp_lat,tmp_lon] = gather(D.time, D.lat, D.lon);

% Generate at least the Number of ADS-B Reports for Normal Operation
ADSB_Spoofing_num = heightT;
ADSB_Spoofing_done = 0;

% ADS-B Spoofing Data
ADSB_Spoofing_idx = zeros(heightT,1,'uint32');
ADSB_Spoofing_icao = zeros(heightT,1,'uint16');
ADSB_Spoofing_lat = zeros(heightT,1);
ADSB_Spoofing_lon = zeros(heightT,1);
ADSB_Spoofing_index_lat = zeros(heightT,1,'uint16');
ADSB_Spoofing_index_lon = zeros(heightT,1,'uint16');
ADSB_Spoofing_sensors = zeros(0,'uint8');   % sensors -> index_sensors?
ADSB_Spoofing_sensors_ref = zeros(heightT,1,'uint16');
ADSB_Spoofing_duration = zeros(heightT,1,'uint16');
ADSB_Spoofing_distance = zeros(heightT,1);

% LUT for Efficient Bit Setting
LUT_ADD = zeros(numel(Sensors),size(index_sensors,2),'uint8');
for i = 1:size(index_sensors,2)
    for j = 1:8
        LUT_ADD((i-1)*8+j,i) = bitset(LUT_ADD((i-1)*8+j,i),j);
        
        % Catch End
        if (i-1)*8+j == numel(Sensors)
            break
        end
    end
end

while ADSB_Spoofing_done < ADSB_Spoofing_num
    tic
    
    % Random Aircraft Report
    select = randi(heightT,1);
    
    % Indices for Reports within the Next Hour
    ADSB_Spoofing_idx_tmp = find(index_icao == index_icao(select) & ...
                                 tmp_time >= tmp_time(select) & ...
                                 tmp_time <= (tmp_time(select) + 3600));
    
    % Duration
    ADSB_Spoofing_duration(ADSB_Spoofing_done+1:ADSB_Spoofing_done+nnz(ADSB_Spoofing_idx_tmp)) = ...
        tmp_time(ADSB_Spoofing_idx_tmp) - tmp_time(ADSB_Spoofing_idx_tmp(1));
    
    % Distances
    ADSB_Spoofing_distance(ADSB_Spoofing_done+1:ADSB_Spoofing_done+nnz(ADSB_Spoofing_idx_tmp)) = ...
        distance(tmp_lat(ADSB_Spoofing_idx_tmp(1)), tmp_lon(ADSB_Spoofing_idx_tmp(1)), ...
                 tmp_lat(ADSB_Spoofing_idx_tmp), tmp_lon(ADSB_Spoofing_idx_tmp), E);
    
    % ADS-B Spoofing Reports
    ADSB_Spoofing_idx(ADSB_Spoofing_done+1:ADSB_Spoofing_done+nnz(ADSB_Spoofing_idx_tmp)) = ADSB_Spoofing_idx_tmp;
    ADSB_Spoofing_icao(ADSB_Spoofing_done+1:ADSB_Spoofing_done+nnz(ADSB_Spoofing_idx_tmp)) = index_icao(ADSB_Spoofing_idx_tmp);
    ADSB_Spoofing_lat(ADSB_Spoofing_done+1:ADSB_Spoofing_done+nnz(ADSB_Spoofing_idx_tmp)) = tmp_lat(ADSB_Spoofing_idx_tmp);
    ADSB_Spoofing_lon(ADSB_Spoofing_done+1:ADSB_Spoofing_done+nnz(ADSB_Spoofing_idx_tmp)) = tmp_lon(ADSB_Spoofing_idx_tmp);
    ADSB_Spoofing_index_lat(ADSB_Spoofing_done+1:ADSB_Spoofing_done+nnz(ADSB_Spoofing_idx_tmp)) = index_lat(ADSB_Spoofing_idx_tmp);
    ADSB_Spoofing_index_lon(ADSB_Spoofing_done+1:ADSB_Spoofing_done+nnz(ADSB_Spoofing_idx_tmp)) = index_lon(ADSB_Spoofing_idx_tmp);
    
    % Sensors that cover the Grid Index
    Sensors_covered = find(Grid(ADSB_Spoofing_index_lat(ADSB_Spoofing_done+1)).Sensors(:,ADSB_Spoofing_index_lon(ADSB_Spoofing_done+1)));
    
    % Mean number of sensors that receive a report within Grid Index
    Sensors_selected = randsample(Sensors_covered, ...
        round(sum(Grid(ADSB_Spoofing_index_lat(ADSB_Spoofing_done+1)).Sensors(:,ADSB_Spoofing_index_lon(ADSB_Spoofing_done+1))) / ...
              double(Grid(ADSB_Spoofing_index_lat(ADSB_Spoofing_done+1)).Total_Messages(:,ADSB_Spoofing_index_lon(ADSB_Spoofing_done+1)))));
    
    % Single
    ADSB_Spoofing_sensors(end+1,:,1) = sum(LUT_ADD(Sensors_selected(1),:),1,'native');                                                
    % Multiple
    ADSB_Spoofing_sensors(end,:,2) = sum(LUT_ADD(Sensors_selected,:),1,'native');
    % All
    ADSB_Spoofing_sensors(end,:,3) = sum(LUT_ADD(Sensors_covered,:),1,'native');
    
    % Ref to Sensors Covered
    ADSB_Spoofing_sensors_ref(ADSB_Spoofing_done+1:ADSB_Spoofing_done+nnz(ADSB_Spoofing_idx_tmp)) = ...
        repmat(size(ADSB_Spoofing_sensors,1),nnz(ADSB_Spoofing_idx_tmp),1);
    
    % Reports Generated
    ADSB_Spoofing_done = ADSB_Spoofing_done + nnz(ADSB_Spoofing_idx_tmp);
    
    fprintf("Successfully generated ADS-B Spoofing reports! %d out of %d (%.2f%%) [Time: %s]\n", ...
            ADSB_Spoofing_done, heightT, 100 * ADSB_Spoofing_done / heightT, seconds(toc));
end
clear ADSB_Spoofing_idx_tmp

% Unique Grid Indices
[grid_index_unique,~,idx_tmp] = unique([ADSB_Spoofing_index_lat, ADSB_Spoofing_index_lon], 'rows');

ADSB_Spoofing_done = 0;
% Save for each Unique Index
for i = 1:size(grid_index_unique,1)
    tic
    
    tmp = find(idx_tmp==i);
    % Table for Spoofing Reports
    ADSB_Spoofing_Single = table(tmp_time(ADSB_Spoofing_idx(tmp)), ...
                                 ADSB_Spoofing_icao(tmp), ...
                                 ADSB_Spoofing_lat(tmp), ...
                                 ADSB_Spoofing_lon(tmp), ...
                                 ADSB_Spoofing_sensors(ADSB_Spoofing_sensors_ref(tmp),:,1), ...
                                 ADSB_Spoofing_duration(tmp), ...
                                 ADSB_Spoofing_distance(tmp), ...
                                 'VariableNames', {'time','icao24','lat','lon','sensors','duration','distance'});
    
    ADSB_Spoofing_Multiple = table(tmp_time(ADSB_Spoofing_idx(tmp)), ...
                                   ADSB_Spoofing_icao(tmp), ...
                                   ADSB_Spoofing_lat(tmp), ...
                                   ADSB_Spoofing_lon(tmp), ...
                                   ADSB_Spoofing_sensors(ADSB_Spoofing_sensors_ref(tmp),:,2), ...
                                   ADSB_Spoofing_duration(tmp), ...
                                   ADSB_Spoofing_distance(tmp), ...
                                   'VariableNames', {'time','icao24','lat','lon','sensors','duration','distance'});
    
    ADSB_Spoofing_All = table(tmp_time(ADSB_Spoofing_idx(tmp)), ...
                              ADSB_Spoofing_icao(tmp), ...
                              ADSB_Spoofing_lat(tmp), ...
                              ADSB_Spoofing_lon(tmp), ...
                              ADSB_Spoofing_sensors(ADSB_Spoofing_sensors_ref(tmp),:,3), ...
                              ADSB_Spoofing_duration(tmp), ...
                              ADSB_Spoofing_distance(tmp), ...
                              'VariableNames', {'time','icao24','lat','lon','sensors','duration','distance'});
    
    % Save Data
    writetable(ADSB_Spoofing_Single, fullfile(filedir_Grid, num2str(grid_index_unique(i,1)), num2str(grid_index_unique(i,2)), "ADSB_Spoofing_Single.csv"));
    writetable(ADSB_Spoofing_Multiple, fullfile(filedir_Grid, num2str(grid_index_unique(i,1)), num2str(grid_index_unique(i,2)), "ADSB_Spoofing_Multiple.csv"));
    writetable(ADSB_Spoofing_All, fullfile(filedir_Grid, num2str(grid_index_unique(i,1)), num2str(grid_index_unique(i,2)), "ADSB_Spoofing_All.csv"));
    
    % Reports Saved
    ADSB_Spoofing_done = ADSB_Spoofing_done + numel(tmp);
    
    fprintf("Successfully saved ADS-B Spoofing reports! %d out of %d (%.2f%%) [Time: %s | Remaining: %s]\n", ...
            ADSB_Spoofing_done, size(ADSB_Spoofing_idx,1), 100 * ADSB_Spoofing_done / size(ADSB_Spoofing_idx,1), ...
            seconds(toc), duration(0,0,toc*(size(grid_index_unique,1)-i)));
end

%%%%%%%%%%%%%%%%%%%%
%%%%% TRAINING %%%%%
%%%%%%%%%%%%%%%%%%%%
fprintf("\nTraining of ADS-B Reports for each Grid Index ...\n");

% LUT for Efficient Reverting of Sensor Indices
LUT_REV = zeros(256,8,'logical');
for i = 1:size(LUT_REV,1)
    for j = 1:8
        LUT_REV(i,j) = bitget(i-1,j);
    end
end

% Number of Sensor Index Bytes
num_idx = ceil(numel(Sensors)/8);

% Load Training Code
addpath(fullfile("..", day));

% For all Grid Indices
for lat = 1:numel(grid_lat)
    tic
    
    for lon = 1:numel(grid_lon{lat})
        % Filepath
        filedir = fullfile(filedir_Grid, num2str(lat), num2str(lon));
        
        % Skip if Already Trained
        %if exist(fullfile(filedir, "trained_Normal_GPS_Fine_Tree.mat"), "file") || ...
        %   exist(fullfile(filedir, "trained_Normal_ADSB_Single_Fine_Tree.mat"), "file")
        %    continue
        %end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% NORMAL OPERATION %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if exist(fullfile(filedir, "Normal_Operation.csv"), "file")
            % Read in Normal Operation Reports
            Normal_Operation = readtable(fullfile(filedir, "Normal_Operation.csv"), ...
                                         'Format',"%u %u16 %f %f " + join(repmat("%u8",1,num_idx)));
            
            % Recover Sensor Index -> Logical Array
            tmp_sensors = table2array(Normal_Operation(:,5:end));
            Normal_Operation_Sensors = false(height(Normal_Operation),numel(Sensors));
            for i = 1:height(Normal_Operation)
                for j = 1:floor(numel(Sensors)/8)
                    Normal_Operation_Sensors(i,(j-1)*8+1:j*8) = LUT_REV(uint16(tmp_sensors(i,j))+1,:);
                end
                % If Number of Sensors is not a Multiple of 8
                if floor(numel(Sensors)/8) ~= num_idx
                    Normal_Operation_Sensors(i,(floor(numel(Sensors)/8)*8)+1:end) = ...
                        LUT_REV(tmp_sensors(i,num_idx)+1,1:(numel(Sensors)-floor(numel(Sensors)/8)*8));
                end
            end
        else
            continue
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% GPS SPOOFING %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%
        if exist(fullfile(filedir, "GPS_Spoofing.csv"), "file")
            % Read in GPS Spoofing Reports
            GPS_Spoofing = readtable(fullfile(filedir, "GPS_Spoofing.csv"), ...
                                     'Format',"%u %u16 %f %f " + join(repmat("%u8",1,num_idx)) + " %u8 %u16 %f");
            
            % Limit Number of Reports to Normal
            if height(GPS_Spoofing) > height(Normal_Operation)
                GPS_Spoofing = GPS_Spoofing(randperm(height(GPS_Spoofing),height(Normal_Operation)),:);
            end
            
            % Recover Sensor Index -> Logical Array
            tmp_sensors = table2array(GPS_Spoofing(:,5:end-3));
            GPS_Spoofing_Sensors = false(height(GPS_Spoofing),numel(Sensors));
            for i = 1:height(GPS_Spoofing)
                for j = 1:floor(numel(Sensors)/8)
                    GPS_Spoofing_Sensors(i,(j-1)*8+1:j*8) = LUT_REV(uint16(tmp_sensors(i,j))+1,:);
                end
                % If Number of Sensors is not a Multiple of 8
                if floor(numel(Sensors)/8) ~= num_idx
                    GPS_Spoofing_Sensors(i,(floor(numel(Sensors)/8)*8)+1:end) = ...
                        LUT_REV(tmp_sensors(i,num_idx)+1,1:(numel(Sensors)-floor(numel(Sensors)/8)*8));
                end
            end
            
            % Perform Training
            Training_Normal_GPS = [table(Normal_Operation_Sensors, ...
                                         repmat(categorical("Normal Operation"),size(Normal_Operation_Sensors,1),1), ...
                                         'VariableNames', {'sensors','Status'}); ...
                                   table(GPS_Spoofing_Sensors, ...
                                         repmat(categorical("GPS Spoofing"),size(GPS_Spoofing_Sensors,1),1), ...
                                         'VariableNames', {'sensors','Status'})];
            trained_Normal_GPS_Fine_Tree = train_Normal_GPS_Fine_Tree(Training_Normal_GPS);
            %trained_Normal_GPS_SVM_Linear = train_Normal_GPS_SVM_Linear(Training_Normal_GPS);
            %trained_Normal_GPS_Fine_KNN = train_Normal_GPS_Fine_KNN(Training_Normal_GPS);
            %trained_Normal_GPS_Boosted_Trees = train_Normal_GPS_Boosted_Trees(Training_Normal_GPS);
            
            % Save Model
            save(fullfile(filedir, "trained_Normal_GPS_Fine_Tree"),"trained_Normal_GPS_Fine_Tree", '-v7.3');
            %save(fullfile(filedir, "trained_Normal_GPS_SVM_Linear"),"trained_Normal_GPS_SVM_Linear", '-v7.3');
            %save(fullfile(filedir, "trained_Normal_GPS_Fine_KNN"),"trained_Normal_GPS_Fine_KNN", '-v7.3');
            %save(fullfile(filedir, "trained_Normal_GPS_Boosted_Trees"),"trained_Normal_GPS_Boosted_Trees", '-v7.3');
            
            % Free Memory
            clear Training_Normal_GPS
            clear trained_Normal_GPS_Fine_Tree
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% ADS-B SPOOFING - SINGLE %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if exist(fullfile(filedir, "ADSB_Spoofing_Single.csv"), "file")
            % Read in ADSB Spoofing Reports
            ADSB_Spoofing_Single = readtable(fullfile(filedir, "ADSB_Spoofing_Single.csv"), ...
                                             'Format',"%u %u16 %f %f " + join(repmat("%u8",1,num_idx)) + " %u16 %f");
            
            % Limit Number of Reports to Normal
            if height(ADSB_Spoofing_Single) > height(Normal_Operation)
                ADSB_Spoofing_Single = ADSB_Spoofing_Single(randperm(height(ADSB_Spoofing_Single),height(Normal_Operation)),:);
            end
            
            % Recover Sensor Index -> Logical Array
            tmp_sensors = table2array(ADSB_Spoofing_Single(:,5:end-2));
            ADSB_Spoofing_Single_Sensors = false(height(ADSB_Spoofing_Single),numel(Sensors));
            for i = 1:height(ADSB_Spoofing_Single)
                for j = 1:floor(numel(Sensors)/8)
                    ADSB_Spoofing_Single_Sensors(i,(j-1)*8+1:j*8) = LUT_REV(uint16(tmp_sensors(i,j))+1,:);
                end
                % If Number of Sensors is not a Multiple of 8
                if floor(numel(Sensors)/8) ~= num_idx
                    ADSB_Spoofing_Single_Sensors(i,(floor(numel(Sensors)/8)*8)+1:end) = ...
                        LUT_REV(tmp_sensors(i,num_idx)+1,1:(numel(Sensors)-floor(numel(Sensors)/8)*8));
                end
            end
            
            % Perform Training
            Training_Normal_ADSB_Single = [table(Normal_Operation_Sensors, ...
                                                 repmat(categorical("Normal Operation"),size(Normal_Operation_Sensors,1),1), ...
                                                 'VariableNames', {'sensors','Status'}); ...
                                           table(ADSB_Spoofing_Single_Sensors, ...
                                                 repmat(categorical("ADSB Spoofing"),size(ADSB_Spoofing_Single_Sensors,1),1), ...
                                                 'VariableNames', {'sensors','Status'})];
            trained_Normal_ADSB_Single_Fine_Tree = train_Normal_ADSB_Fine_Tree(Training_Normal_ADSB_Single);
            %trained_Normal_ADSB_Single_SVM_Linear = train_Normal_ADSB_SVM_Linear(Training_Normal_ADSB_Single);
            %trained_Normal_ADSB_Single_Fine_KNN = train_Normal_ADSB_Fine_KNN(Training_Normal_ADSB_Single);
            %trained_Normal_ADSB_Single_Boosted_Trees = train_Normal_ADSB_Boosted_Trees(Training_Normal_ADSB_Single);
            
            % Save Model
            save(fullfile(filedir, "trained_Normal_ADSB_Single_Fine_Tree"),"trained_Normal_ADSB_Single_Fine_Tree", '-v7.3');
            %save(fullfile(filedir, "trained_Normal_ADSB_Single_SVM_Linear"),"trained_Normal_ADSB_Single_SVM_Linear", '-v7.3');
            %save(fullfile(filedir, "trained_Normal_ADSB_Single_Fine_KNN"),"trained_Normal_ADSB_Single_Fine_KNN", '-v7.3');
            %save(fullfile(filedir, "trained_Normal_ADSB_Single_Boosted_Trees"),"trained_Normal_ADSB_Single_Boosted_Trees", '-v7.3');
            
            % Free Memory
            clear Training_Normal_ADSB_Single
            clear trained_Normal_ADSB_Single_Fine_Tree
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% ADS-B SPOOFING - MULTIPLE %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if exist(fullfile(filedir, "ADSB_Spoofing_Multiple.csv"), "file")
            % Read in ADSB Spoofing Reports
            ADSB_Spoofing_Multiple = readtable(fullfile(filedir, "ADSB_Spoofing_Multiple.csv"), ...
                                               'Format',"%u %u16 %f %f " + join(repmat("%u8",1,num_idx)) + " %u16 %f");
            
            % Limit Number of Reports to Normal
            if height(ADSB_Spoofing_Multiple) > height(Normal_Operation)
                ADSB_Spoofing_Multiple = ADSB_Spoofing_Multiple(randperm(height(ADSB_Spoofing_Multiple),height(Normal_Operation)),:);
            end
                                           
            % Recover Sensor Index -> Logical Array
            tmp_sensors = table2array(ADSB_Spoofing_Multiple(:,5:end-2));
            ADSB_Spoofing_Multiple_Sensors = false(height(ADSB_Spoofing_Multiple),numel(Sensors));
            for i = 1:height(ADSB_Spoofing_Multiple)
                for j = 1:floor(numel(Sensors)/8)
                    ADSB_Spoofing_Multiple_Sensors(i,(j-1)*8+1:j*8) = LUT_REV(uint16(tmp_sensors(i,j))+1,:);
                end
                % If Number of Sensors is not a Multiple of 8
                if floor(numel(Sensors)/8) ~= num_idx
                    ADSB_Spoofing_Multiple_Sensors(i,(floor(numel(Sensors)/8)*8)+1:end) = ...
                        LUT_REV(tmp_sensors(i,num_idx)+1,1:(numel(Sensors)-floor(numel(Sensors)/8)*8));
                end
            end
            
            % Perform Training
            Training_Normal_ADSB_Multiple = [table(Normal_Operation_Sensors, ...
                                                   repmat(categorical("Normal Operation"),size(Normal_Operation_Sensors,1),1), ...
                                                   'VariableNames', {'sensors','Status'}); ...
                                             table(ADSB_Spoofing_Multiple_Sensors, ...
                                                   repmat(categorical("ADSB Spoofing"),size(ADSB_Spoofing_Multiple_Sensors,1),1), ...
                                                   'VariableNames', {'sensors','Status'})];
            trained_Normal_ADSB_Multiple_Fine_Tree = train_Normal_ADSB_Fine_Tree(Training_Normal_ADSB_Multiple);
            %trained_Normal_ADSB_Multiple_SVM_Linear = train_Normal_ADSB_SVM_Linear(Training_Normal_ADSB_Multiple);
            %trained_Normal_ADSB_Multiple_Fine_KNN = train_Normal_ADSB_Fine_KNN(Training_Normal_ADSB_Multiple);
            %trained_Normal_ADSB_Multiple_Boosted_Trees = train_Normal_ADSB_Boosted_Trees(Training_Normal_ADSB_Multiple);
            
            % Save Model
            save(fullfile(filedir, "trained_Normal_ADSB_Multiple_Fine_Tree"),"trained_Normal_ADSB_Multiple_Fine_Tree", '-v7.3');
            %save(fullfile(filedir, "trained_Normal_ADSB_Multiple_SVM_Linear"),"trained_Normal_ADSB_Multiple_SVM_Linear", '-v7.3');
            %save(fullfile(filedir, "trained_Normal_ADSB_Multiple_Fine_KNN"),"trained_Normal_ADSB_Multiple_Fine_KNN", '-v7.3');
            %save(fullfile(filedir, "trained_Normal_ADSB_Multiple_Boosted_Trees"),"trained_Normal_ADSB_Multiple_Boosted_Trees", '-v7.3');
            
            % Free Memory
            clear Training_Normal_ADSB_Multiple
            clear trained_Normal_ADSB_Multiple_Fine_Tree
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% ADS-B SPOOFING - ALL %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if exist(fullfile(filedir, "ADSB_Spoofing_All.csv"), "file")
            % Read in ADSB Spoofing Reports
            ADSB_Spoofing_All = readtable(fullfile(filedir, "ADSB_Spoofing_All.csv"), ...
                                          'Format',"%u %u16 %f %f " + join(repmat("%u8",1,num_idx)) + " %u16 %f");
            
            % Limit Number of Reports to Normal
            if height(ADSB_Spoofing_All) > height(Normal_Operation)
                ADSB_Spoofing_All = ADSB_Spoofing_All(randperm(height(ADSB_Spoofing_All),height(Normal_Operation)),:);
            end
            
            % Recover Sensor Index -> Logical Array
            tmp_sensors = table2array(ADSB_Spoofing_All(:,5:end-2));
            ADSB_Spoofing_All_Sensors = false(height(ADSB_Spoofing_All),numel(Sensors));
            for i = 1:height(ADSB_Spoofing_All)
                for j = 1:floor(numel(Sensors)/8)
                    ADSB_Spoofing_All_Sensors(i,(j-1)*8+1:j*8) = LUT_REV(uint16(tmp_sensors(i,j))+1,:);
                end
                % If Number of Sensors is not a Multiple of 8
                if floor(numel(Sensors)/8) ~= num_idx
                    ADSB_Spoofing_All_Sensors(i,(floor(numel(Sensors)/8)*8)+1:end) = ...
                        LUT_REV(tmp_sensors(i,num_idx)+1,1:(numel(Sensors)-floor(numel(Sensors)/8)*8));
                end
            end
            
            % Perform Training
            Training_Normal_ADSB_All = [table(Normal_Operation_Sensors, ...
                                              repmat(categorical("Normal Operation"),size(Normal_Operation_Sensors,1),1), ...
                                              'VariableNames', {'sensors','Status'}); ...
                                        table(ADSB_Spoofing_All_Sensors, ...
                                              repmat(categorical("ADSB Spoofing"),size(ADSB_Spoofing_All_Sensors,1),1), ...
                                              'VariableNames', {'sensors','Status'})];
            trained_Normal_ADSB_All_Fine_Tree = train_Normal_ADSB_Fine_Tree(Training_Normal_ADSB_All);
            %trained_Normal_ADSB_All_SVM_Linear = train_Normal_ADSB_SVM_Linear(Training_Normal_ADSB_All);
            %trained_Normal_ADSB_All_Fine_KNN = train_Normal_ADSB_Fine_KNN(Training_Normal_ADSB_All);
            %trained_Normal_ADSB_All_Boosted_Trees = train_Normal_ADSB_Boosted_Trees(Training_Normal_ADSB_All);
            
            % Save Model
            save(fullfile(filedir, "trained_Normal_ADSB_All_Fine_Tree"),"trained_Normal_ADSB_All_Fine_Tree", '-v7.3');
            %save(fullfile(filedir, "trained_Normal_ADSB_All_SVM_Linear"),"trained_Normal_ADSB_All_SVM_Linear", '-v7.3');
            %save(fullfile(filedir, "trained_Normal_ADSB_All_Fine_KNN"),"trained_Normal_ADSB_All_Fine_KNN", '-v7.3');
            %save(fullfile(filedir, "trained_Normal_ADSB_All_Boosted_Trees"),"trained_Normal_ADSB_All_Boosted_Trees", '-v7.3');
            
            % Free Memory
            clear Training_Normal_ADSB_All
            clear trained_Normal_ADSB_All_Fine_Tree
        end
        
        %%%%%%%%%%%%%%%%%%
        %%%%% HYBRID %%%%%
        %%%%%%%%%%%%%%%%%%
        if exist(fullfile(filedir, "GPS_Spoofing.csv"), "file") && ...
           exist(fullfile(filedir, "ADSB_Spoofing_Multiple.csv"), "file")
            % Perform Training
            Training_Normal_GPS_ADSB = [table(Normal_Operation_Sensors, ...
                                              repmat(categorical("Normal Operation"),size(Normal_Operation_Sensors,1),1), ...
                                              'VariableNames', {'sensors','Status'}); ...
                                        table(GPS_Spoofing_Sensors, ...
                                              repmat(categorical("GPS Spoofing"),size(GPS_Spoofing_Sensors,1),1), ...
                                              'VariableNames', {'sensors','Status'}); ...
                                        table(ADSB_Spoofing_Multiple_Sensors, ...
                                              repmat(categorical("ADSB Spoofing"),size(ADSB_Spoofing_Multiple_Sensors,1),1), ...
                                              'VariableNames', {'sensors','Status'})];
            trained_Normal_GPS_ADSB_Fine_Tree = train_Normal_GPS_ADSB_Fine_Tree(Training_Normal_GPS_ADSB);
            %trained_Normal_GPS_ADSB_SVM_Linear = train_Normal_GPS_ADSB_SVM_Linear(Training_Normal_GPS_ADSB);
            %trained_Normal_GPS_ADSB_Fine_KNN = train_Normal_GPS_ADSB_Fine_KNN(Training_Normal_GPS_ADSB);
            %trained_Normal_GPS_ADSB_Boosted_Trees = train_Normal_GPS_ADSB_Boosted_Trees(Training_Normal_GPS_ADSB);
            
            % Save Model
            save(fullfile(filedir, "trained_Normal_GPS_ADSB_Fine_Tree"),"trained_Normal_GPS_ADSB_Fine_Tree", '-v7.3');
            %save(fullfile(filedir, "trained_Normal_GPS_ADSB_SVM_Linear"),"trained_Normal_GPS_ADSB_SVM_Linear", '-v7.3');
            %save(fullfile(filedir, "trained_Normal_GPS_ADSB_Fine_KNN"),"trained_Normal_GPS_ADSB_Fine_KNN", '-v7.3');
            %save(fullfile(filedir, "trained_Normal_GPS_ADSB_Boosted_Trees"),"trained_Normal_GPS_ADSB_Boosted_Trees", '-v7.3');
            
            % Free Memory
            clear Training_Normal_GPS_ADSB
            clear trained_Normal_GPS_ADSB_Fine_Tree
        end
    end
    fprintf("Successful Training! Latitude Index %d out of %d (%.2f%%) [Time: %s | Remaining: %s]\n", ...
            lat, numel(grid_lat), 100 * lat / numel(grid_lat), ...
            seconds(toc), duration(0,0,toc*(numel(grid_lat)-lat)));
end
fprintf("Training successfully finished!\n");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% TESTING DATA GENERATION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf("\nTesting Data Generation: Importing Reports ...\n");

% Datastore
ds = tabularTextDatastore(fullfile(filedir_Reports, "*.csv"));
ds.SelectedFormats = report_format;
ds.SelectedVariableNames = {'time','lat','lon'};
D = tall(ds);

% Read Time Latitude Longitude
[tmp_time,tmp_lat,tmp_lon] = gather(D.time, D.lat, D.lon);

% 1000 Runs
run_num = 1000;

% Scenario (1) Normal Operation
%          (2) GPS Spoofing - Deviation 1°
%          (3) GPS Spoofing - Deviation 2°
%          (4) GPS Spoofing - Deviation 5°
%          (5) GPS Spoofing - Deviation 10°
%          (6) GPS Spoofing - Deviation 20°
%          (7) GPS Spoofing - Deviation 45°
%          (8) ADSB Spoofing - Single
%          (9) ADSB Spoofing - Multiple
%         (10) ADSB Spoofing - All
%         (11) GPS + ADSB Spoofing - Deviation 5° + Multiple

% GPS Spoofing Deviation
deviation = [1, 2, 5, 10, 20, 45];

% Timw Window
t_before = 900; % 15 Mins before
t_after = 3600; % 60 Mins after

% Testing Data [run x scenario x time]
Data_Testing_index_lat = zeros(run_num,11,t_before+1+t_after,'uint16');
Data_Testing_index_lon = zeros(run_num,11,t_before+1+t_after,'uint16');
Data_Testing_index_sensors = zeros(run_num,11,t_before+1+t_after,ceil(numel(Sensors)/8),'uint8');
Data_Testing_distance = zeros(run_num,11,t_before+1+t_after);

% Sample Tracks
run_done = 0;
while run_done < run_num
    tic
    
    % Random Aircraft Report
    select = randi(heightT,1);
    
    % Check if Track 15 Mins before and 60 Mins after exists
    if isempty(find(index_icao == index_icao(select) & tmp_time == (tmp_time(select) - t_before), 1)) || ...
       isempty(find(index_icao == index_icao(select) & tmp_time == (tmp_time(select) + t_after), 1))
        continue
    end    
    
    % Indices for Reports
    Data_idx_tmp = find(index_icao == index_icao(select) & ...
                        tmp_time >= (tmp_time(select) - t_before)& ...
                        tmp_time <= (tmp_time(select) + t_after));
    
    % Check if at least 3600 Reports
    if numel(Data_idx_tmp) < t_after
        continue
    end
    
    % Reports for Track
    track = table(tmp_time(Data_idx_tmp), ...
                  index_icao(Data_idx_tmp), ...
                  tmp_lat(Data_idx_tmp), ...
                  tmp_lon(Data_idx_tmp), ...
                  index_sensors(Data_idx_tmp,:), ...
                  'VariableNames', {'time','icao24','lat','lon','sensors'});
    
    % Duration
    Data_duration_tmp = track.time - track.time(1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% NORMAL OPERATION %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    Data_Testing_index_lat(run_done+1,1,Data_duration_tmp+1) = index_lat(Data_idx_tmp);
    Data_Testing_index_lon(run_done+1,1,Data_duration_tmp+1) = index_lon(Data_idx_tmp);
    
    Data_Testing_index_sensors(run_done+1,1,Data_duration_tmp+1,:) = index_sensors(Data_idx_tmp,:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% GPS SPOOFING %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%
    % Distances and Angles to all the next Points (when attack starts)
    [arclen,az] = distance(track.lat(track.time==tmp_time(select)), track.lon(track.time==tmp_time(select)), ...
                           track.lat(track.time>=tmp_time(select)), track.lon(track.time>=tmp_time(select)), E);
    
    % For all Deviations
    for dev = 1:numel(deviation)
        
        % Calculate new Positions with Deviation of x° (only plus)
        [GPS_Spoofing_lat_tmp,GPS_Spoofing_lon_tmp] = reckon(track.lat(track.time==tmp_time(select)), track.lon(track.time==tmp_time(select)), ...
                                                             arclen, az+deviation(dev), E);
        
        % Calculate Grid Index
        GPS_Spoofing_index_lat_tmp = zeros(numel(GPS_Spoofing_lat_tmp),1,'uint16');
        GPS_Spoofing_index_lon_tmp = zeros(numel(GPS_Spoofing_lon_tmp),1,'uint16');
        for i = 1:numel(GPS_Spoofing_lat_tmp)
            % Distance to all Lat in Grid
            [~,GPS_Spoofing_index_lat_tmp(i)] = min(distance('rh', GPS_Spoofing_lat_tmp(i), 0, ...
                                                                   grid_lat, 0, E));
            % Distance to all Lon in Grid
            [~,GPS_Spoofing_index_lon_tmp(i)] = min(distance('rh', Grid(GPS_Spoofing_index_lat_tmp(i)).lat, GPS_Spoofing_lon_tmp(i), ...
                                                                   Grid(GPS_Spoofing_index_lat_tmp(i)).lat, Grid(GPS_Spoofing_index_lat_tmp(i)).lon, E));
        end
        
        % GPS Spoofing Reports
        Data_Testing_index_lat(run_done+1,dev+1,1:t_before) = Data_Testing_index_lat(run_done+1,1,1:t_before);
        Data_Testing_index_lat(run_done+1,dev+1,Data_duration_tmp(Data_duration_tmp>=t_before)+1) = GPS_Spoofing_index_lat_tmp;
        
        Data_Testing_index_lon(run_done+1,dev+1,1:t_before) = Data_Testing_index_lon(run_done+1,1,1:t_before);
        Data_Testing_index_lon(run_done+1,dev+1,Data_duration_tmp(Data_duration_tmp>=t_before)+1) = GPS_Spoofing_index_lon_tmp;
        
        Data_Testing_index_sensors(run_done+1,dev+1,:,:) = Data_Testing_index_sensors(run_done+1,1,:,:);
        
        Data_Testing_distance(run_done+1,dev+1,Data_duration_tmp(Data_duration_tmp>=t_before)+1) = ...
            distance(track.lat(track.time>=tmp_time(select)), track.lon(track.time>=tmp_time(select)), ...
                     GPS_Spoofing_lat_tmp, GPS_Spoofing_lon_tmp, E);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% ADS-B SPOOFING %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % LUT for Efficient Bit Setting
    LUT_ADD = zeros(numel(Sensors),size(index_sensors,2),'uint8');
    for i = 1:size(index_sensors,2)
        for j = 1:8
            LUT_ADD((i-1)*8+j,i) = bitset(LUT_ADD((i-1)*8+j,i),j);
            
            % Catch End
            if (i-1)*8+j == numel(Sensors)
                break
            end
        end
    end  
    
    % Sensors that cover the Grid Index
    Sensors_covered = find(Grid(index_lat(Data_idx_tmp(Data_duration_tmp==t_before))).Sensors(:,index_lon(Data_idx_tmp(Data_duration_tmp==t_before))));
    
    % Mean number of sensors that receive a report within Grid Index
    Sensors_selected = randsample(Sensors_covered, ...
        round(sum(Grid(index_lat(Data_idx_tmp(Data_duration_tmp==t_before))).Sensors(:,index_lon(Data_idx_tmp(Data_duration_tmp==t_before)))) / ...
              double(Grid(index_lat(Data_idx_tmp(Data_duration_tmp==t_before))).Total_Messages(:,index_lon(Data_idx_tmp(Data_duration_tmp==t_before))))));
    
    % Single
    ADSB_Spoofing_sensors(1,:) = sum(LUT_ADD(Sensors_selected(1),:),1,'native');
    % Multiple
    ADSB_Spoofing_sensors(2,:) = sum(LUT_ADD(Sensors_selected,:),1,'native');
    % All
    ADSB_Spoofing_sensors(3,:) = sum(LUT_ADD(Sensors_covered,:),1,'native');
    
    % ADS-B Spoofing Reports    
    Data_Testing_index_lat(run_done+1,8:10,:) = repmat(Data_Testing_index_lat(run_done+1,1,:),1,3);
    Data_Testing_index_lon(run_done+1,8:10,:) = repmat(Data_Testing_index_lon(run_done+1,1,:),1,3);
    
    Data_Testing_index_sensors(run_done+1,8:10,1:t_before,:) = repmat(Data_Testing_index_sensors(run_done+1,1,1:t_before,:),1,3);
    Data_Testing_index_sensors(run_done+1,8:10,Data_duration_tmp(Data_duration_tmp>=t_before)+1,:) = ...
        permute(repmat(ADSB_Spoofing_sensors,1,1,nnz(Data_duration_tmp(Data_duration_tmp>=t_before))),[1,3,2]);
    
    Data_Testing_distance(run_done+1,8:10,Data_duration_tmp(Data_duration_tmp>=t_before)+1) = ...
        repmat(distance(track.lat(Data_duration_tmp==t_before), track.lon(Data_duration_tmp==t_before), ...
                        track.lat(track.time>=tmp_time(select)), track.lon(track.time>=tmp_time(select)), E),1,3)';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% GPS + ADS-B SPOOFING %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate new Positions with Deviation of 5° (only plus)
    [GPS_Spoofing_lat_tmp,GPS_Spoofing_lon_tmp] = reckon(track.lat(track.time==tmp_time(select)), track.lon(track.time==tmp_time(select)), ...
                                                         arclen, az+5, E);
    
    % Calculate Grid Index
    GPS_Spoofing_index_lat_tmp = zeros(numel(GPS_Spoofing_lat_tmp),1,'uint16');
    GPS_Spoofing_index_lon_tmp = zeros(numel(GPS_Spoofing_lon_tmp),1,'uint16');
    for i = 1:numel(GPS_Spoofing_lat_tmp)
        % Distance to all Lat in Grid
        [~,GPS_Spoofing_index_lat_tmp(i)] = min(distance('rh', GPS_Spoofing_lat_tmp(i), 0, ...
                                                               grid_lat, 0, E));
        % Distance to all Lon in Grid
        [~,GPS_Spoofing_index_lon_tmp(i)] = min(distance('rh', Grid(GPS_Spoofing_index_lat_tmp(i)).lat, GPS_Spoofing_lon_tmp(i), ...
                                                               Grid(GPS_Spoofing_index_lat_tmp(i)).lat, Grid(GPS_Spoofing_index_lat_tmp(i)).lon, E));
    end
    
    % GPS + ADS-B Spoofing Reports
    Data_Testing_index_lat(run_done+1,11,1:t_before) = Data_Testing_index_lat(run_done+1,1,1:t_before);
    Data_Testing_index_lat(run_done+1,11,Data_duration_tmp(Data_duration_tmp>=t_before)+1) = GPS_Spoofing_index_lat_tmp;
    
    Data_Testing_index_lon(run_done+1,11,1:t_before) = Data_Testing_index_lon(run_done+1,1,1:t_before);
    Data_Testing_index_lon(run_done+1,11,Data_duration_tmp(Data_duration_tmp>=t_before)+1) = GPS_Spoofing_index_lon_tmp;
    
    Data_Testing_index_sensors(run_done+1,11,:,:) = Data_Testing_index_sensors(run_done+1,9,:,:);
    
    Data_Testing_distance(run_done+1,11,Data_duration_tmp(Data_duration_tmp>=t_before)+1) = ...
        distance(track.lat(track.time==tmp_time(select)), track.lon(track.time==tmp_time(select)), ...
                 GPS_Spoofing_lat_tmp, GPS_Spoofing_lon_tmp, E);
    
    % Track Generated
    run_done = run_done + 1;
    
    fprintf("Successfully generated reports for testing! %d out of %d (%.2f%%) [Time: %s | Remaining: %s]\n", ...
            run_done, run_num, 100 * run_done / run_num, seconds(toc), duration(0,0,toc*(run_num-run_done)));
end

% Save Testing Data
save(fullfile(filedir_Data, "Data_Testing.mat"), "Data_Testing_index_lat", "Data_Testing_index_lon", "Data_Testing_index_sensors", "Data_Testing_distance", '-v7.3')

% Free Memory
clear Data_idx_tmp
clear Data_duration_tmp
clear GPS_Spoofing_lat_tmp
clear GPS_Spoofing_lon_tmp
clear GPS_Spoofing_index_lat_tmp
clear GPS_Spoofing_index_lon_tmp
clear index_icao
clear index_lat
clear index_lon
clear index_sensors

%%%%%%%%%%%%%%%%%%%
%%%%% TESTING %%%%%
%%%%%%%%%%%%%%%%%%%
fprintf("\nTesting of ADS-B Reports for each Grid Index ...\n");

% Load Testing Data
load(fullfile(filedir_Data, "Data_Testing.mat"))

% Results for GPS | ADSB_Single | ADSB_Multiple | ADSB_All | GPS + ADSB_Mulitple
%Results = repmat(categorical("Not Classified"),1000,11*4501,5);
load(fullfile(filedir_Data, "Results.mat"))

% LUT for Efficient Reverting of Sensor Indices
LUT_REV = zeros(256,8,'logical');
for i = 1:size(LUT_REV,1)
    for j = 1:8
        LUT_REV(i,j) = bitget(i-1,j);
    end
end

% Number of Sensor Index Bytes
num_idx = ceil(numel(Sensors)/8);

% For all Runs
for run = 1:1000
    tic
    
    % Reshape Data Testing
    Data_Testing_index_lat_tmp = reshape(Data_Testing_index_lat(run,:,:),11*4501,1);
    Data_Testing_index_lon_tmp = reshape(Data_Testing_index_lon(run,:,:),11*4501,1);
    Data_Testing_index_sensors_tmp = reshape(Data_Testing_index_sensors(run,:,:,:),11*4501,num_idx);
    
    % Unique Grid Indices
    [grid_index_unique,~,idx_tmp] = unique([Data_Testing_index_lat_tmp, Data_Testing_index_lon_tmp], 'rows');
    
    % For each Unique Index
    for index = 1:size(grid_index_unique,1)
        
        % Skip [0|0]
        if grid_index_unique(index,1) == 0 && grid_index_unique(index,2) == 0
            continue
        end
        
        % Filepath
        filedir = fullfile(filedir_Grid, num2str(grid_index_unique(index,1)), num2str(grid_index_unique(index,2)));
        
        % Recover Sensor Index -> Logical Array
        tmp_sensors = squeeze(Data_Testing_index_sensors_tmp(idx_tmp==index,:));         
        
        Data_Sensors = false(size(tmp_sensors,1),numel(Sensors));
        for i = 1:size(tmp_sensors,1)
            for j = 1:floor(numel(Sensors)/8)
                Data_Sensors(i,(j-1)*8+1:j*8) = LUT_REV(uint16(tmp_sensors(i,j))+1,:);
            end
            % If Number of Sensors is not a Multiple of 8
            if floor(numel(Sensors)/8) ~= num_idx
                Data_Sensors(i,(floor(numel(Sensors)/8)*8)+1:end) = ...
                    LUT_REV(tmp_sensors(i,num_idx)+1,1:(numel(Sensors)-floor(numel(Sensors)/8)*8));
            end
        end
        
        % Testing Data
        Testing = table(Data_Sensors, 'VariableNames', {'sensors'});
        
        % Check if Model (Normal|GPS) for Grid Index Exists
        if exist(fullfile(filedir, "trained_Normal_GPS_Fine_Tree.mat"), "file")
            % Read in Model
            load(fullfile(filedir, "trained_Normal_GPS_Fine_Tree.mat"))
            
            % Testing Normal + GPS
            Results(run,idx_tmp==index,1) = trained_Normal_GPS_Fine_Tree.predictFcn(Testing);
            
            % Free Memory
            clear trained_Normal_GPS_Fine_Tree
        end
        
        % Check if Model (Normal|ADSB_Single) for Grid Index Exists
        if exist(fullfile(filedir, "trained_Normal_ADSB_Single_Fine_Tree.mat"), "file")
            % Read in Model
            load(fullfile(filedir, "trained_Normal_ADSB_Single_Fine_Tree.mat"))
            
            % Testing Normal + ADSB
            Results(run,idx_tmp==index,2) = trained_Normal_ADSB_Single_Fine_Tree.predictFcn(Testing);
            
            % Free Memory
            clear trained_Normal_ADSB_Single_Fine_Tree
        end
        
        % Check if Model (Normal|ADSB_Multiple) for Grid Index Exists
        if exist(fullfile(filedir, "trained_Normal_ADSB_Multiple_Fine_Tree.mat"), "file")
            % Read in Model
            load(fullfile(filedir, "trained_Normal_ADSB_Multiple_Fine_Tree.mat"))
            
            % Testing Normal + ADSB
            Results(run,idx_tmp==index,3) = trained_Normal_ADSB_Multiple_Fine_Tree.predictFcn(Testing);
            
            % Free Memory
            clear trained_Normal_ADSB_Multiple_Fine_Tree
        end
        
        % Check if Model (Normal|ADSB_All) for Grid Index Exists
        if exist(fullfile(filedir, "trained_Normal_ADSB_All_Fine_Tree.mat"), "file")
            % Read in Model
            load(fullfile(filedir, "trained_Normal_ADSB_All_Fine_Tree.mat"))
            
            % Testing Normal + ADSB
            Results(run,idx_tmp==index,4) = trained_Normal_ADSB_All_Fine_Tree.predictFcn(Testing);
            
            % Free Memory
            clear trained_Normal_ADSB_All_Fine_Tree
        end
        
        % Check if Model (Normal|GPS|ADSB_Multiple) for Grid Index Exists
        if exist(fullfile(filedir, "trained_Normal_GPS_ADSB_Fine_Tree.mat"), "file")
            % Read in Model
            load(fullfile(filedir, "trained_Normal_GPS_ADSB_Fine_Tree.mat"))
            
            % Testing Normal + GPS + ADSB
            Results(run,idx_tmp==index,5) = trained_Normal_GPS_ADSB_Fine_Tree.predictFcn(Testing);
            
            % Free Memory
            clear trained_Normal_GPS_ADSB_Fine_Tree
        end
    end     
    fprintf("Successful Testing! %d out of %d (%.2f%%) [Time: %s | Remaining: %s]\n", ...
            run, 1000, 100 * run / 1000, ...
            seconds(toc), duration(0,0,toc*(1000-run)));
end

% Save Results
save(fullfile(filedir_Data, "Results.mat"), "Results", '-v7.3')

%%%%%%%%%%%%%%%%%%%%%
%%%%% EVALUATON %%%%%
%%%%%%%%%%%%%%%%%%%%%
fprintf("\nEvaluation of Attack Detection Performance ...\n");

% Load Results
load(fullfile(filedir_Data, "Results.mat"))
load(fullfile(filedir_Data, "Data_Testing.mat"))

run_num = 1000;

% Reshape Results [RUN x SCENARIO x TIME x CLASSIFIER]
Results_reshaped = reshape(Results(1:run_num,:,:),run_num,11,4501,5);

% Convert Results
Results_converted = double(Results_reshaped);
%%% MAY CHANGE - CHECK %%%
% (1) Not Classified
% (2) GPS Spoofing
% (3) Normal Operation
% (4) ADSB Spoofing

% Set "Not Classified" to Normal Operation or NaN? (Currently to Normal)
Results_converted(Results_converted==1) = 3;
Results_converted(Results_converted==1) = NaN;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% GPS SPOOFING - TABLE - TIME WINDOWS - THRESHOLDS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Everything Normal (Even for Spoofing)
idx = ~ismember(squeeze(Results_converted(:,2,:,1)), ones(1,4501)*3, 'rows');

% ALPHA [1,2,5,10,20,45]
% TIME WINDOWS 10, 15, 20 Mins
Results_smoothed = NaN(run_num,10,4501,3);
w = [300, 600, 900];
for i = 1:numel(w)
    for scenario = 1:11
        % Smooth Results with moving average
        Results_smoothed(:,scenario,:,i) = movmean(squeeze(Results_converted(:,scenario,:,1)),[w(i) 0],2);
    end
end

% Threshold (Minimal Score for "Normal Operation")
threshold = min(squeeze(min(Results_smoothed(1:run_num,1,901:end,:))));
%threshold_half = min(squeeze(min(Results_smoothed(1:run_num,1,901:2701,:))));

% False Alarm Rate
for alpha = 1:6
    for w = 1:3
        alert = [];
        for run = 1:run_num
            if ~idx(run)
                continue
            end

            alert = [alert, find(Results_smoothed(run,alpha+1,901:end,w) < threshold(w),1) - 1];
        end

        fprintf("ALPHA: " + alpha + " - W: " + w + " - Alerts: " + numel(alert) + " (" + numel(alert)/nnz(idx) + ...
            ") - T_MEAN: " + mean(alert)/60 + " - T_MEDIAN: " + median(alert)/60 + " - T_STD: " + std(alert/60) + "\n");
    end
end

%Data_Testing_distance(11,7,2701)

dist_res = 0:100:30000;
res = NaN(300,1);
for i = 2:numel(dist_res)
    idx = Data_Testing_distance(1:run_num,2:7,:)>dist_res(i-1) & Data_Testing_distance(1:run_num,2:7,:)<=dist_res(i);

    tmp = squeeze(Results_converted(:,2:7,:,1));
    
    res(i-1) = mean(tmp(idx));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% GRID COVERAGE FOR TEST DATA %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data_Testing_Grid_Coverage = NaN(size(Data_Testing_distance));
for run = 1:171
    for scenario = 1:10
        for t = 1:4501
            if ~Data_Testing_index_lat(run,scenario,t)
                continue
            end
            
            Data_Testing_Grid_Coverage(run,scenario,t) = nnz(Grid(Data_Testing_index_lat(run,scenario,t)).Sensors(:,Data_Testing_index_lon(run,scenario,t)));
        end
    end
end


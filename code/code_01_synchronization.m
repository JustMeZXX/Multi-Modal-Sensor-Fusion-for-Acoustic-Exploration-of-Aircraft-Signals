clc;clear;close all

load('labels_2021217_2021247.mat') % load labels

ini_Jday = 216; % initial Jday
ini_Mday = 4; % day in August for Jday 216

day_a_month = 31;

% parameters
c = 330; % m/s
reflat = 32.84591;
reflon = -96.78411;
refalt = 131;

for file_ind = 1:day_a_month
    
    clc;close all;
    clearvars -except gt_segment_all ini_Jday ini_Mday day_a_month...
        c reflat reflon refalt file_ind
    
    %% date setup
    cur_day = ini_Mday + file_ind;
    
    if cur_day <= day_a_month
        time_ini = sprintf('2021-08-%02d 00:00:00',cur_day);
    else
        time_ini = sprintf('2021-09-%02d 00:00:00',mod(cur_day,day_a_month));
        cur_day = mod(cur_day,day_a_month);
    end
    
    t_ini = datevec(datenum(time_ini));
    time_interval = 60*60*24;
    
    %% load .csv file
    cur_Jday = ini_Jday + file_ind;
    csv_name = sprintf('adsb_2021%03d.csv',cur_Jday);
    
    adsb_table_all = readtable(csv_name);
    
    %% filter out duplicate rows with seen_pos > 1s
    seen_pos = 1;
    ind = find(table2array(adsb_table_all(:,14)) <= seen_pos);
    adsb_table = adsb_table_all(ind,:);
    
    %% filter out NaN rows
    id_nan = isnan(table2array(adsb_table(:,8))) |...
        isnan(table2array(adsb_table(:,9))) |...
        isnan(table2array(adsb_table(:,10))) |...
        isnan(table2array(adsb_table(:,11))) |...
        isnan(table2array(adsb_table(:,12))) |...
        isnan(table2array(adsb_table(:,28)));
    
    adsb_table = adsb_table(~id_nan,:);
    
    %% find valid index of time interval in adsb file before time shifting
    time_for_pick_num = datenum(datetime(table2array(adsb_table(:,4)),"ConvertFrom","epochtime"));
    [~, sort_idx] = sort(time_for_pick_num);
    
    adsb_table = adsb_table(sort_idx,:);
    
    time_for_pick_date = datetime(table2array(adsb_table(:,4)),"ConvertFrom","epochtime");
    [~,~,day] = ymd(time_for_pick_date);
    
    ind_start = find(day == cur_day,1);
    
    if isempty(ind_start)
        ind_start = 1;
    end
    
    ind_end = find(day == (cur_day+1),1);
    
    if isempty(ind_end)
        ind_end = size(adsb_table,1);
    end
    
    adsb_table = adsb_table(ind_start:ind_end,:);
    
    %% cluster aircrafts following aircraft hex
    [aircraft_hex,~] = unique(adsb_table(:,5));
    
    adsb_table_cluster = [];
    for kplane = 1:size(aircraft_hex,1)
        
        kkpind = find(strcmp(aircraft_hex{kplane,1},adsb_table{:,5}));
        adsb_table_cluster{kplane} = adsb_table(kkpind,:);
    end
    
    %% synchronization
    gt_segment = gt_segment_all{file_ind}; % day 217 to 247 (index 1 to 31)
    gt_to_adsb_table = [];
    
    for gt_ind = 1:size(gt_segment,1)
        
        gt_segment_cur = gt_segment(gt_ind,1:2);
        for k = 1:length(adsb_table_cluster)
            
            segment_positive_for_cur_gt_table = [];
            time_sequence_cur = zeros(1,time_interval);
            
            adsb_table_cluster_cur = adsb_table_cluster{k};
            
            %% calculate time difference
            alt = adsb_table_cluster_cur{:,8} * 0.3048; % feet to m
            lat = adsb_table_cluster_cur{:,9};
            lon = adsb_table_cluster_cur{:,10};
            track = adsb_table_cluster_cur{:,11};
            vel = adsb_table_cluster_cur{:,12} * 0.44704; % miles/h to m/s
            
            [arclen,az] = distance(lat,lon,reflat*ones(size(lat)),reflon*ones(size(lon)),'degrees');
            
            dist = deg2km(arclen) * 1000; % m
            dist3 = sqrt(dist.^2 + (alt-refalt).^2);
            tilt = atan((alt-refalt)./dist);
            
            vdot = vel .* cos((track-az)/180*pi) .* cos(tilt);
            t_delta = dist3 ./ (c + vdot); % time shift for an aircraft event
            
            %% compensate the time shift (t_delta)
            adsb_table_cluster_cur.message_date = adsb_table_cluster_cur{:,4};
            
            t_cur = datevec(datetime(table2array(adsb_table_cluster_cur(:,4))+t_delta,"ConvertFrom","epochtime")); % ignore column 14 due to rounding accuracy
            time_diff = floor(etime(t_cur,t_ini));
            t_val_ind = find(time_diff >= 0 & time_diff < time_interval-1);
            
            if isempty(t_val_ind)
                gt_to_adsb_table{gt_ind}{k} = [];
                continue
            end
            
            adsb_table_cluster_cur{t_val_ind,3} = table2array(adsb_table_cluster_cur(t_val_ind,4))+t_delta(t_val_ind);
            adsb_table_cluster_cur = adsb_table_cluster_cur(t_val_ind,:);
            
            %% process the synchronized cluster
            time_for_pick_num = datenum(datetime(table2array(adsb_table_cluster_cur(:,3)),"ConvertFrom","epochtime"));
            [~, sort_idx] = sort(time_for_pick_num);
            
            adsb_table_cluster_cur = adsb_table_cluster_cur(sort_idx,:);
            
            t_cur = datevec(datetime(table2array(adsb_table_cluster_cur(:,3)),"ConvertFrom","epochtime"));
            time_diff = floor(etime(t_cur,t_ini));
            
            time_sequence_cur(time_diff+1) = time_sequence_cur(time_diff+1) + 1;
            time_diff_cur = time_diff + 1; % map of time diff
            
            time_sequence_cur_binary = time_sequence_cur > 0;
            
            time_sequence_cur_binary_in_gt = time_sequence_cur_binary(1,gt_segment_cur(1):gt_segment_cur(2));
            time_sequence_cur_binary_in_gt_tmp = [0,time_sequence_cur_binary_in_gt]; % add dummy 0
            
            %% find adsb segments with binary 1 in gt
            segment_positive = [];
            positive_index = find(diff(time_sequence_cur_binary_in_gt_tmp')==1);
            
            for i = 1:length(positive_index)
                
                index_start_cur = positive_index(i);
                index_start_cur_temp = index_start_cur;
                
                while index_start_cur_temp <= length(time_sequence_cur_binary_in_gt) && time_sequence_cur_binary_in_gt(index_start_cur_temp) == 1
                    index_start_cur_temp = index_start_cur_temp + 1;
                end
                
                index_end_cur = index_start_cur_temp - 1;
                index_cur = [index_start_cur,index_end_cur];
                
                segment_positive = [segment_positive;int32(index_cur)]; % returned segments from adsb
            end
            
            if ~isempty(segment_positive)
                segment_positive_final = segment_positive + gt_segment_cur(1) - 1;
            else
                segment_positive_final = [];
            end
            
            for kk = 1:size(segment_positive_final,1)
                
                [~, start_ind] = min(abs(time_diff_cur - double(segment_positive_final(kk,1))));
                [~, end_ind] = min(abs(time_diff_cur - double(segment_positive_final(kk,2))));
                
                segment_positive_for_cur_gt_table = [segment_positive_for_cur_gt_table; adsb_table_cluster_cur(start_ind:end_ind,:)];
            end
            
            if ~isempty(segment_positive_for_cur_gt_table)
                
                [~,unique_ind] = unique(segment_positive_for_cur_gt_table(:,1));
                
                segment_positive_for_cur_gt_table_uni = segment_positive_for_cur_gt_table(unique_ind,:);
            else
                segment_positive_for_cur_gt_table_uni = [];
            end
            
            gt_to_adsb_table{gt_ind}{k} = segment_positive_for_cur_gt_table_uni;
        end
    end
    
    % after synchronization
    % the result (gt_to_adsb_table) is organized as a 2-D cell
    % the first index iterates the corresponding manual labels
    % the second index iterates the aircraft events
    
    % for example,
    % for the ith label, the jth aircraft event is gt_to_adsb_table{i}{j}
    
    gt_to_adsb_result.table = gt_to_adsb_table;
    
    save_name = sprintf('gt_to_adsb_syn_result_%03d.mat',cur_Jday);
    save(save_name,'gt_to_adsb_result')
end


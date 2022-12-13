clc;clear;close all

load('labels_2021217_2021247.mat') % load labels

ini_Jday = 216; % initial Jday
ini_Mday = 4; % day in August for Jday 216

day_a_month = 31;

% columns with useful aircraft information
valid_cols = 1:49;
save_cols = [3,4,5,7,8,9,10,11,12,14,28,49,50,51];

% parameters
c = 330;
reflat = 32.84591;
reflon = -96.78411;
refalt = 131;

% pre-defined thresholds
range_threshold = 32000; % unit in m
doppler_shift_scale_threshold = 0.2; % ratio
tk_margin = 10; % unit in s

for file_ind = 1:day_a_month
    
    clc;close all;
    clearvars -except gt_segment_all ini_Jday ini_Mday day_a_month...
        valid_cols save_cols c reflat reflon refalt...
        range_threshold doppler_shift_scale_threshold tk_margin file_ind
    
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
    
    %% load synchronized results
    cur_Jday = ini_Jday + file_ind;
    mat_name = sprintf('gt_to_adsb_syn_result_%03d.mat',cur_Jday);
    
    load(mat_name);
    
    gt_segment = gt_segment_all{file_ind};
    gt_to_adsb_table_all = gt_to_adsb_result.table;
    
    %% filter out aircraft not contributing to Doppler shift
    gt_to_adsb_table = [];
    
    for iigt = 1:size(gt_to_adsb_table_all,2)
        for iiair = 1:size(gt_to_adsb_table_all{iigt},2)
            
            gt_to_adsb_table_temp = gt_to_adsb_table_all{iigt}{iiair};
            if ~isempty(gt_to_adsb_table_temp)
                
                lat = gt_to_adsb_table_temp{:,9};
                lon = gt_to_adsb_table_temp{:,10};
                track = gt_to_adsb_table_temp{:,11};
                
                vel = gt_to_adsb_table_temp{:,12} * 0.44704; % m/h to m/s
                
                [arclen,az] = distance(lat,lon,reflat*ones(size(lat)),reflon*ones(size(lon)),'degrees');
                
                dist = deg2km(arclen) * 1000;
                alt = gt_to_adsb_table_temp{:,8} * 0.3048;
                dist3 = sqrt(dist.^2 + (alt-refalt).^2);
                
                tilt = atan((alt-refalt)./dist);
                
                vdot = vel .* cos((track-az)/180*pi) .* cos(tilt);
                
                doppler_shifit = c ./ (c - vdot);
                
                ia = find(dist3 < range_threshold);
                iaa = ia(~isoutlier(doppler_shifit(ia)));
                
                if isempty(iaa) || (max(doppler_shifit(iaa)) - min(doppler_shifit(iaa))) < doppler_shift_scale_threshold
                    gt_to_adsb_table{iigt}{iiair} = [];
                else
                    gt_to_adsb_table{iigt}{iiair} = gt_to_adsb_table_temp(iaa,:);
                end
            else
                gt_to_adsb_table{iigt}{iiair} = [];
            end
        end
    end
    
    %% filter out aircraft not aligning with tk (appearing time)
    for iigt = 1:size(gt_to_adsb_table,2)
        
        tk_cur = gt_segment(iigt,4);
        for iiair = 1:size(gt_to_adsb_table{iigt},2)
            
            gt_to_adsb_table_cur = gt_to_adsb_table{iigt}{iiair};
            if ~isempty(gt_to_adsb_table_cur)
                
                t_cur = datevec(datetime(table2array(gt_to_adsb_table_cur(:,3)),"ConvertFrom","epochtime"));
                time_diff = floor(etime(t_cur,t_ini))+1;
                
                % find whether time of aircraft at closest point in t_cur
                % is shifted with tk within a margin (10s)
                
                lat = gt_to_adsb_table_cur{:,9};
                lon = gt_to_adsb_table_cur{:,10};
                
                [arclen,az] = distance(lat,lon,reflat*ones(size(lat)),reflon*ones(size(lon)),'degrees');
                
                dist = deg2km(arclen) * 1000; % m
                alt = gt_to_adsb_table_cur{:,8} * 0.3048; % feet to m
                dist3 = sqrt(dist.^2 + (alt-refalt).^2);
                
                [dist3_min, dist3_min_index] = min(dist3);
                
                diff_tk_cur = abs(tk_cur - time_diff(dist3_min_index));
                
                if (diff_tk_cur > tk_margin)
                    gt_to_adsb_table{iigt}{iiair} = [];
                else
                    
                    % append the maximum distance of aircraft
                    gt_to_adsb_table_cur.('maximum_distance') = max(dist3) * ones(size(gt_to_adsb_table_cur,1),1);
                    gt_to_adsb_table{iigt}{iiair} = gt_to_adsb_table_cur;
                end
            end
        end
    end
    
    adsb_into_types_gt = [];
    
    for gt_ind = 1:length(gt_to_adsb_table)
        
        gt_segment_cur = gt_segment(gt_ind,1:2);
        
        gt_to_aircraft_cur = zeros(length(gt_to_adsb_table{gt_ind}),time_interval);
        gt_to_type_cur = zeros(length(gt_to_adsb_table{gt_ind}),time_interval);
        
        count = 0;
        recored_index = [];
        
        for aircraft_ind = 1:length(gt_to_adsb_table{gt_ind})
            
            if ~isempty(gt_to_adsb_table{gt_ind}{aircraft_ind})
                
                count = count + 1;
                
                gt_to_adsb_table_cur = gt_to_adsb_table{gt_ind}{aircraft_ind};
                
                time_sequence_cur = zeros(1,time_interval);
                time_diff_cur = [];
                
                for t = 1:size(gt_to_adsb_table_cur,1)
                    
                    t_cur = datevec(datetime(table2array(gt_to_adsb_table_cur(t,3)),"ConvertFrom","epochtime"));
                    time_diff = floor(etime(t_cur,t_ini));
                    
                    if (time_diff >= 0) && (time_diff < time_interval-1)
                        time_sequence_cur(time_diff+1) = time_sequence_cur(time_diff+1) + 1;
                        
                        time_diff_cur = [time_diff_cur;time_diff + 1]; % map of time diff
                    end
                end
                
                % record time index of adsb to each aircraft regarding gt
                recored_index{aircraft_ind} = time_diff_cur;
                
                gt_to_aircraft_cur(count,time_diff_cur) = aircraft_ind;
                if ~isnan(unique(gt_to_adsb_table{gt_ind}{aircraft_ind}{:,28}))
                    
                    uni_tmp = unique(gt_to_adsb_table{gt_ind}{aircraft_ind}{:,28});
                    gt_to_type_cur(count,time_diff_cur) = uni_tmp; % aircraft type
                    
                else
                    gt_to_type_cur(count,time_diff_cur) = 0; % 0 to replace NaN
                end
            end
        end
        
        %% Find index of columns with different aircraft types (0 is presaved)
        gt_to_type_cur_seg = gt_to_type_cur(:,gt_segment_cur(1):gt_segment_cur(2));
        
        base_types = {4,5}; % FWSE and FWME
        
        % predefine the target types for visualization
        base_types_ind_gt_cur = [];
        
        for iii = 1:length(base_types)
            
            base_types_cur = base_types{iii};
            base_types_ind_cur = [];
            
            for jjj = 1:size(gt_to_type_cur_seg,2)
                
                gt_to_type_cur_seg_temp = gt_to_type_cur_seg(:,jjj);
                gt_to_type_cur_seg_temp(gt_to_type_cur_seg_temp==0) = [];
                gt_to_type_cur_seg_temp_uni = unique(gt_to_type_cur_seg_temp);
                
                if isequal(gt_to_type_cur_seg_temp_uni,base_types_cur)
                    
                    % conver jjj to global index
                    base_types_ind_cur = [base_types_ind_cur,jjj+gt_segment_cur(1)-1];
                end
            end
            
            base_types_ind_gt_cur{iii} = base_types_ind_cur;
        end
        
        % merge the adsb table into different types (only for format)
        adsb_into_types = [];
        
        for types_ind = 1:length(base_types_ind_gt_cur)
            
            types_ind_cur = base_types_ind_gt_cur{types_ind};
            gt_to_aircraft_ind = gt_to_aircraft_cur(:,types_ind_cur);
            
            gt_to_aircraft_ind_all = [];
            adsb_into_types_temp = [];
            
            for iid = 1:size(gt_to_aircraft_ind,2)
                
                gt_to_aircraft_ind_cur = unique(gt_to_aircraft_ind(:,iid));
                gt_to_aircraft_ind_cur = gt_to_aircraft_ind_cur(gt_to_aircraft_ind_cur~=0);
                
                gt_to_aircraft_ind_all = [gt_to_aircraft_ind_all; gt_to_aircraft_ind_cur];
            end
            
            gt_to_aircraft_ind_uni = unique(gt_to_aircraft_ind_all);
            for iij = 1:length(gt_to_aircraft_ind_uni)
                
                % append source frequency and its appearing time
                freq_cur = gt_segment(gt_ind, 3);
                time_cur = gt_segment(gt_ind, 4);
                
                for iik = 1:length(types_ind_cur)
                    
                    indkk = find(recored_index{gt_to_aircraft_ind_uni(iij)} == types_ind_cur(iik));
                    
                    table_cur = gt_to_adsb_table{gt_ind}{gt_to_aircraft_ind_uni(iij)}(indkk,valid_cols);
                    table_cur.('source_freq') = double(freq_cur) * ones(size(table_cur,1),1);
                    table_cur.('appearing_time') = double(time_cur) * ones(size(table_cur,1),1);
                    
                    adsb_into_types_temp = [adsb_into_types_temp; table_cur];
                end
            end
            
            adsb_into_types{types_ind} = adsb_into_types_temp;
        end
        
        adsb_into_types_gt{gt_ind} = adsb_into_types;
    end
    
    % after association
    % the result (adsb_into_types_day_cur) is organized as a cell of table
    % the index iterates the aircraft type, (1 for FWSE and 2 for FWME)
    % the resultant adsb data for each type is stored as table in the cell
    
    % for example,
    % the adsb data (in table) for FWSE is adsb_into_types_day_cur{1}
    % the adsb data (in table) for FWME is adsb_into_types_day_cur{2}
    
    adsb_into_types_day_temp = [];
    for kkk = 1:size(adsb_into_types_gt,2)
        adsb_into_types_day_temp = [adsb_into_types_day_temp;adsb_into_types_gt{kkk}];
    end
    
    adsb_into_types_day_cur = [];
    for lll = 1:size(adsb_into_types_day_temp,2)
        
        adsb_into_types_day_cur_temp = cat(1,adsb_into_types_day_temp{:,lll});
        
        if ~isempty(adsb_into_types_day_cur_temp)
            adsb_into_types_day_cur{lll} = adsb_into_types_day_cur_temp(:,save_cols);
        else
            adsb_into_types_day_cur{lll} = [];
        end
    end
    
    save_name = sprintf('gt_to_adsb_assoc_result_%05d_%03d_%03d.mat',range_threshold,tk_margin,cur_Jday);
    save(save_name,'adsb_into_types_day_cur')
end


%% THREE POINT BENDING PROCESSING PIPELINE
% Written for Blaine Christiansen's Musculoskeletal Adaptation Lab
% University of California Davis Health
% Department of Orthopaedic Surgery
% By Dovin Kiernan
% 2025 Feb 12

%% PURPOSE
% Takes .txt files output by ???THREE POINT BENDING MACHINE??? ???PROTOCOL???
% Batch process them using methods based on Jepsen et al. 2015
% Allow user to review and edit the automatic processing
% Output an Excel file with the following variables:
%%% 1) Filename
%%% 2) Stiffness (N/mm)
%%% 3) yield displacement (mm)
%%% 4) yield load (N)
%%% 5) post-yield displacement (mm)
%%% 6) fracture displacement (mm)
%%% 7) fracture load (N)
%%% 8) max load (N)
%%% 9) work to yield (N*mm)
%%% 10) work to fracture (N*mm)
%%% 11) whether trial was manually edited (0 or 1)
%%% 12) whether trial was rejected (0 or 1)

%% CLEAN WORKSPACE
clear; close; clc;

%% SET CONSTANTS
% Size between frames in mm
% Smaller is noisier but more resolution for finding the yield point, I've been using 0.001 or 0.0001 mm
interp_size = 0.0001;
% If any break points are being found too early and being erroneously
% defined as the end of the linear portion of the data,
% you can increase the threshold to decrease the number of break points found
linear_threshold = 2;
% Displacement changes ~ between -2e-3 and + 2e-3 then has a very sharp change when fracture occurs
% In the example data this was > 0.4 and was very clear, so I used a threshold of 0.1
fracture_threshold = 0.1;
% UCD colours
poppy = [241/256 138/256 0];
strawberry = [249/256 53/256 73/256];
putahcreek = [0 142/256 170/256];
% Declare global variables to share with call back (button press) functions
global stiffness stiffness_line yield_idx yield_line replot reject_trial manual_trial

%% UI OPEN FOLDER CONTAINING FILES
data_folder = uigetdir;
cd(data_folder)
% Populate a list of files to be processed
files_for_processing = dir;
% Remove everything except 'asc' files
for file_count = size(files_for_processing,1):-1:1
    if files_for_processing(file_count).isdir || ~strcmp(files_for_processing(file_count).name(end-2:end),'TXT')
        files_for_processing(file_count) = [];
    end
end % file_count
% Make table to store output
output = table;
% Create folder to store processed data
mkdir([data_folder,'/Processed'])

%% LOAD AND LOOP THROUGH EACH FILE
for file_count = 1:size(files_for_processing,1)
    reject_trial = 0; manual_trial = 0;
    % Specify import options
    opts = detectImportOptions([data_folder,'\',files_for_processing(file_count).name],...
        'delimiter',',', ...
        'DecimalSeparator','.', ...
        'VariableNamesLine',4);
    % Import
    data_raw = readtable([data_folder,'\',files_for_processing(file_count).name],opts);
    % Some of the example trials were empty, which would crash the code, so skip and reject empty trials
    if isempty(data_raw)
        % Write outputs with filename
        output.name{file_count} = files_for_processing(file_count).name;
        output.stiffness(file_count) = NaN; % Correct from N/sample to N/mm
        output.yield_displacement(file_count) = NaN; % mm
        output.yield_load(file_count) = NaN; % N
        output.postyield_displacement(file_count) = NaN; % mm
        output.fracture_displacement(file_count) = NaN; % mm
        output.fracture_load(file_count) = NaN; % N
        output.max_load(file_count) = NaN; % N
        output.work_to_yield(file_count) = NaN;
        output.work_to_fracture(file_count) = NaN;
        output.manual(file_count) = 0;
        output.rejected(file_count) = 1;
    else
        % Clear first row (this is where the units were written)
        data_raw(1,:) = [];
        % Should now have...
        %%% "Points" -- this is sample *within* a 1 s scan and goes from 1 to 200 then resets, don't need, remove
        data_raw.Points = [];
        %%% "ElapsedTime" -- timestamp in s (0.005 s per example data provided by Sophie Orr), continuous throughout collection with no resetting
        %%% every 200 samples a new scan occurs and an empty row appears in the data, remove
        for row_count = 201:200:size(data_raw,1)-floor(size(data_raw,1)/200)
            data_raw(row_count,:) = [];
        end
        %%% "ScanTime" -- timestamp in s (0.005 s per sample default), resets every 1 s scan, don't need, remove
        data_raw.ScanTime = [];
        %%% "Disp1" -- displacement in mm, multiply by -1 to make displacement positive
        data_raw.Disp1 = data_raw.Disp1.*-1;
        %%% "Load1" -- force in N, multiply by -1 to make force positive
        data_raw.Load1 = data_raw.Load1.*-1;
        %%% "Load2" -- not used, remove
        data_raw.Load2 = [];
        % Identify section of file where displacement is occuring
        % Start of trial?
        % Example trials are all starting from point where displacement begins increasing
        % I assume that this the protocol begins recording and applying the displacement simultaneously
        % Therefore, I haven't tried to auto-identify trial start here
        % End of trial --> find point of fracture
        % Displacement changes ~ between -2e-3 and + 2e-3 then has a very sharp change when fracture occurs
        % Find this sharp change in displacement
        [fracture_displacement, fracture_load] = find(diff(data_raw.Disp1)>fracture_threshold,1,'first');
        % Remove data after fracture
        data_raw(fracture_displacement-1:end,:) = [];
        % Data should now no longer have a reversal in displacement values and they should ~basically increase monotonically (in displacement control)
        % However, the displacement and diff(displacement) are noisy with many repeating (non-unique) values
        % I am guessing this occurs because the rate of displacement and the resolution are lower than the sample frequency so the read out is between values for consecutive samples
        % Therefore, we will interpolate to smooth the data and create all unique values based on the Force-displacement curve
        % To do so, create a vector of random numbers with an extremely low magnitude that won't affect our data
        % e-15 is the smallest value that isn't rounded and still counts as "unique"
        random = ((randperm(size(data_raw,1))-round(size(data_raw,1)/2)).*1e-15)';
        % Then create a table and store interpolated time, displacement, and load
        data_interp = table;
        data_interp.Disp = (min(data_raw.Disp1):interp_size:max(data_raw.Disp1))'; % 0.0001 mm per frame (or whatever "interp_size" is set to)
        % Some load trials have very sharp, high frequency, high magnitude noise
        % I tried Butterworth filters and wavelet filters but found the most effective approach to be
        % Using a 0.25 s moving median find and replace elements >3 local scaled MAD (i.e., a Hampel filter)
        data_interp.Load = filloutliers(interp1(data_raw.Disp1+random,data_raw.Load1,data_interp.Disp,'linear','extrap'),'linear','movmedian',round(0.25/mean(diff(data_raw.ElapsedTime)))); % N
        % data_interp.Load = wdenoise(interp1(data_raw.Disp1+random,data_raw.Load1,data_interp.Disp,'linear','extrap'),floor(log2(size(data_interp.Disp,1))),Wavelet='coif1'); % N
        data_interp.Time = interp1(data_raw.Disp1+random,data_raw.ElapsedTime,data_interp.Disp,'linear','extrap'); % s
        % Find maximum load
        max_load = max(data_interp.Load);
        % Find end of linear portion of curve
        % if any break points are being found too early and being erroneously
        % defined as the end of the linear portion of the data,
        % you can increase the threshold to decrease the number of break points found
        linear_disp_opts = ischange(data_interp.Load,'linear','Threshold',linear_threshold);
        % Example data showed that sometimes stiffness increased throughout the linear region
        % Per Blaine, you take the second "stiffer" part of the curve as the stiffness
        % To execute, look at the slopes of each region until slopes are negative
        % Take the largest slope out of those options
        linear_disp_opts = find(linear_disp_opts == 1);
        linear_disp_opts = [1; linear_disp_opts];
        for region_count = 1:size(linear_disp_opts,1)-1
            slopes(region_count,1) = (data_interp.Load(linear_disp_opts(region_count+1)) - data_interp.Load(linear_disp_opts(region_count)))/...
                (linear_disp_opts(region_count+1) - linear_disp_opts(region_count));
            if slopes(region_count,1) < 0
                slopes(region_count,:) = [];
                break
            end
        end
        % Check to see that the code has found slopes
        if ~isempty(slopes)
            % If it has...
            [~, slope_of_interest] = max(slopes);
            slope_of_interest = [slope_of_interest; slope_of_interest + 1];
            linear_disp_idx = linear_disp_opts(slope_of_interest);
        else
            % Otherwise, choose obviously wrong values
            slopes = 100;
            slope_of_interest = 1;
            linear_disp_idx = [1 2];
        end
        % Find stiffness by calculating slope within selected linear portion
        % Rise over run
        % Units in N/sample (whatever displacement between frames is set at above -- leaving like this for plotting and calculation purposes then changing to true units for output)
        stiffness = slopes(slope_of_interest(1));
        % Find the load when the displacement is zero ("b" in y = mx + b)
        zero_intercept = data_interp.Load(linear_disp_idx(2)) - stiffness*linear_disp_idx(2);
        % Calculate a line starting from the 0-intercept
        stiffness_line = (stiffness*(1:size(data_interp,1)) + zero_intercept)';
        % Subtract 10% from slope and find point of intercept for yield
        yield_slope = stiffness*0.9;
        yield_line = (yield_slope*(1:size(data_interp,1)) + zero_intercept)';
        % Operationalize this as the first point after end_linear_displacement that falls below the yield line
        yield_idx = find((data_interp.Load - yield_line) < 0);
        yield_idx = yield_idx(yield_idx>linear_disp_idx(2));
        yield_idx = yield_idx(1);
        % Plot results
        replot = 1;
        while replot == 1
            processed_parent = figure;
            set(processed_parent,'Units','Normalized','OuterPosition',[0 0 1 1],...
                'Name', 'THREE POINT BENDING PROCESSOR',...
                'ToolBar','none',...
                'Color', [1 1 1]);
            processed_plot = subplot('Position',[0.05 0.1 0.65 0.85]);
            plot(data_interp.Load,'Color',putahcreek); hold on; plot(stiffness_line,'Color',strawberry); plot(yield_line,'Color',poppy); scatter(yield_idx,data_interp.Load(yield_idx),'MarkerEdgeColor',poppy,'MarkerFaceColor',poppy);
            legend({'data','stiffness','yield'},'Location','north','Box','off')
            title(files_for_processing(file_count).name)
            % Accept or reject
            % If rejected, allow user to reselect trial start and end, and choose points to determine stiffness values
            button_accept = uicontrol(processed_parent,...
                'Units','Normalized','Position',[0.75 0.85 0.2 0.1],...
                'Style','pushbutton',...
                'String','ACCEPT',...
                'Callback',{@accept,processed_parent});
            button_manual = uicontrol(processed_parent,...
                'Units','Normalized','Position',[0.75 0.65 0.2 0.1],...
                'Style','pushbutton',...
                'String','EDIT MANUALLY',...
                'Callback',{@manual_edit,data_interp,processed_parent});
            button_reject = uicontrol(processed_parent,...
                'Units','Normalized','Position',[0.75 0.45 0.2 0.1],...
                'Style','pushbutton',...
                'String','REJECT',...
                'Callback',{@reject_fcn,processed_parent});
            uiwait(processed_parent)
        end
        % Calculate variables
        % Find yield load and displacement
        yield_load = data_interp.Load(yield_idx);
        yield_displacment = data_interp.Disp(yield_idx);
        % Find post-yield displacement
        post_yield_displacement = data_interp.Disp(end) - yield_displacment;
        % Find work
        % Integrate the area under the curve
        work = cumtrapz(data_interp.Disp,data_interp.Load);
        % Work to yield
        work_to_yield = work(yield_idx);
        % Work to fracture
        work_to_fracture = work(end);
        % Write outputs with filename
        output.name{file_count} = files_for_processing(file_count).name;
        output.stiffness(file_count) = stiffness/interp_size; % Correct from N/sample to N/mm
        output.yield_displacement(file_count) = yield_displacment; % mm
        output.yield_load(file_count) = yield_load; % N
        output.postyield_displacement(file_count) = post_yield_displacement; % mm
        output.fracture_displacement(file_count) = data_interp.Disp(end); % mm
        output.fracture_load(file_count) = data_interp.Load(end); % N
        output.max_load(file_count) = max_load; % N
        output.work_to_yield(file_count) = work_to_yield;
        output.work_to_fracture(file_count) = work_to_fracture;
        output.rejected(file_count) = manual_trial;
        output.rejected(file_count) = reject_trial;
    end % if empty
end % file_count
cd([data_folder,'/Processed'])
writetable(output,strcat('MATLAB_3PointOutput_',datestr(now,'yy_mm_dd'),'.xlsx'))

%% ACCEPT FUNCTION
% This function will...
% Activate on button press
% Tell the main function to stop plotting the data and proceed
function accept(source,event,processed_parent)
global replot
close(processed_parent)
replot = 0;
end

%% EDIT FUNCTION
% This function will...
% Activate on button press
% Pop up a new window with all data
% Allow user to manually select start and end of trial (ends at fracture)
% Redraw the figure
% Then allow user to manually select the start and end of the linear region
function manual_edit(source,event,data_interp,processed_parent)
global stiffness stiffness_line yield_idx yield_line replot manual_trial
% New figure
reprocessing_figure = figure;
set(reprocessing_figure,'Units','Normalized','OuterPosition',[0.1 0.1 0.8 0.8],...
    'Name', 'MANUAL REPROCESSING',...
    'ToolBar','none',...
    'Color', [1 1 1]);
reprocessing_plot = subplot('Position',[0.1 0.1 0.8 0.8]);
plot(data_interp.Load); hold on
title('CLICK AT START AND END OF LINEAR REGION')
linear_region = ginput(2);
linear_region = round(linear_region);
% Find stiffness by calculating slope within selected linear portion
% Rise over run
% Units in N/sample (whatever displacement between frames is set at above -- leaving like this for plotting and calculation purposes then changing to true units for output)
stiffness = (data_interp.Load(linear_region(2)) - data_interp.Load(linear_region(1)))/...
    (linear_region(2) - linear_region(1));
% Find the load when the displacement is zero ("b" in y = mx + b)
zero_intercept = data_interp.Load(linear_region(2)) - stiffness*linear_region(2);
% Calculate a line starting from the 0-intercept
stiffness_line = (stiffness*(1:size(data_interp,1)) + zero_intercept)';
% Subtract 10% from slope and find point of intercept for yield
yield_slope = stiffness*0.9;
yield_line = (yield_slope*(1:size(data_interp,1)) + zero_intercept)';
% Operationalize this as the first point after end_linear_displacement that falls below the yield line
yield_idx = find((data_interp.Load - yield_line) < 0);
yield_idx = yield_idx(yield_idx>linear_region(2));
yield_idx = yield_idx(1);
% Plot
close(reprocessing_figure)
close(processed_parent)
replot = 1;
manual_trial = 1;
end % function

%% REJECT FUNCTION
% This function will...
% Activate on button press
% Tell the main function to stop plotting the data and proceed
function reject_fcn(source,event,processed_parent)
global replot reject_trial
close(processed_parent)
replot = 0;
reject_trial = 1;
end
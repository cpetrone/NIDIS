% This code is used to extrapolate and export grey values from images. 
% Values are extracted along multiple straight parallel lines. 
% These lines are perpendicular to the so-called guideline which the user 
% must draw along the feature of interest. To trace the guideline, always 
% remember that the first value will be extrapolated 90º to the right 
% looking down the guideline, from the first point drawn.

% The code uses standard deviation (1sigma) to clear outliers. 
% Values outside the range (mean +/- 1sigma) are cleared (LINE 491). 

% Mean values and other basic statistics values (e.g.: min, max, standard 
% deviation, standard error) are also calculated and exported to a CSV file 
% as well as all the greyvalues extrapolated. 

%  Output folder:
%       CSV file with all the greyvalues and relevant calculations
%       CSV file without outliers and relevamt calculation
%       Reference image
%       Workspace with the guideline coordinates.

% The code has been run on machines running Matlab 2011b, 2014b and 2015a and 
% Windows XP, Windows 7, Windows 10, Ubuntu 12.04 and OS 10.7 and 10.10. 
% It requires Image Processing TOOLBOX (Imline and Improfile) and 
% Statistics and Machine Learning TOOLBOX (pdist, nanmean, nanmin, nanmax, nanstd)


function greyvalues()

while true % Allows to perform multiple calculations without running the code again 
    clc; % Next three lines clear matlab working space and close any open images
    clear;
    close all;
     
    % Start calculating the time for the code to run (useful for programming purposes)
    % tic;
    
    % Select image to analyse and create CSV file names to export the
    % values
    [ImageName,ImagePath] = uigetfile({'*.jpg;*.bmp;*.gif;*.png;*.tif','Images (*.jpg,*.bmp,*.gif,*.png,*.tif)'},'Select the IMAGE');
    imshow(strcat(ImagePath, ImageName));
    figure(gcf);
    I=imread(strcat(ImagePath, ImageName));
    
    % Draw or import the guideline (grey values will be
    % extracted along lines perpendicular to the guideline) 
    % NB: if importing guideline coordinates, make sure you're using the same image as the one the coordinates were saved from.  
    
    gline_msg = menu('How would you like to draw the guideline?','Draw manually','Use existing coordinates');
    if gline_msg == 1
        guideline = imline;
        gline_coordinates = wait(guideline);
        
        % Uncomment the following lines to use ginput function instead of imline
        %{
        [x_guideline,y_guideline] = ginput(2);
        l = line(x_guideline,y_guideline);
        gline_coor_msg = 'Is the guideline correct? (y/n)? ';
        gline_coor_msg_input = input(gline_coor_msg, 's'); 
        while gline_coor_msg_input == 'n'
            disp('Select again the start and end points for the guideline');
            [x_guideline,y_guideline] = ginput(2);
            delete(l);
            l = line(x_guideline,y_guideline);
            gline_coor_msg = 'Is the guideline correct? (y/n)? ';
            gline_coor_msg_input = input(gline_coor_msg, 's'); 
        end
        gline_coordinates = horzcat(x_guideline, y_guideline);
        gline_msg = 'Would you like to save the line coordinates (y/n)? ';
        save_coordinates = input(gline_msg, 's'); 
        if save_coordinates == 'y'
            sc_msg = 'Enter a Filename for the line coordinates: ';
            sc_name = input(sc_msg, 's');
            full_var_name = strcat(sc_name, '.mat');
            save(full_var_name, 'gline_coordinates');
        end 
        %}
        
    elseif gline_msg == 2
        [line_name,line_path] = uigetfile('*.mat','Select the file with the line coordinates');
        load(strcat(line_path, line_name));  
        line(gline_coordinates(:,1),gline_coordinates(:,2), 'Marker', 's');
    elseif gline_msg == 0
        return;
    end
    
    % Calculations to trace lines perpendicular to the guideline (i.e.: greyvalue PROFILE lines)
    cat=diff(gline_coordinates); catx=abs(cat(1)); caty=abs(cat(2));
    len = sqrt(catx.^2+caty.^2);
    
    hp = 50; % HALF length of the profile line.
    fp = 2*hp; % FULL length of the profile line. 
    
    % Uncomment the following lines if you want to trace a custom number of profile lines
    % prompt = 'How many profile lines would you like to trace? ';
    % prof_int_num_temp = input(prompt); 
    
    % Trace as many profile lines as the number of pixels in the guideline.
    % Comment the following lines to trace a custom number of profile lines
    prof_int_num_temp = int32(len);
    
    % Calculations to draw profile lines spaced approximately one pixel
    % from the
    prof_int_num = double(prof_int_num_temp) - 1;
    prof_int = len/(prof_int_num);
    offset = 1;
    
    % NB: MATLAB imread coordinates: y increaes downwards, x increases rightwards 
    
    % Confirm the length of the profile line and extrapolate grey values along straight parallel profile lines perpendicular to the guideline.
    % Different 'if' statements are used to ensure the first value extracted is always on the right of the first point drawn when tracing the guideline.
    
    if gline_coordinates(1)<gline_coordinates(2) && gline_coordinates(3)>gline_coordinates(4)
        [pl_1, pl_2, pl_3, pl_4, xi, yi, x0, y0, x1, y1] = lines(gline_coordinates, catx, caty, hp, len); % See function at the end of the file        
        confirm_msg = 'Is the profile line of the right length: ';
        confirm = strcat(confirm_msg, sprintf(' %d pixels (y/n)? ', hp));
        answer_1 = input(confirm, 's');
        while answer_1 == 'n'     
            new_length = 'Enter the half profile length? ';
            hp = input(new_length);
            fp = 2*hp;
            delete(pl_1);
            delete(pl_2);
            delete(pl_3);
            delete(pl_4);
            [pl_1, pl_2, pl_3, pl_4, xi, yi, x0, y0, x1, y1] = lines(gline_coordinates, catx, caty, hp, len);
            confirm_msg_1 = 'Is the profile line of the right length: ';
            confirm_1 = strcat(confirm_msg_1, sprintf(' %d pixels (y/n)? ', hp));
            answer_1 = input(confirm_1, 's'); 
        end
        delete(pl_2);
        delete(pl_3);
        delete(pl_4);        
        profile_all = zeros(fp, prof_int_num_temp); % Create empty matrix for all the profile values
        [cx, cy, profile_3d] = improfile(I,xi,yi,fp);
        hold on;
        plot(x0, y0, 'Marker','o','MarkerFaceColor','red');
        hold off;
        Vprofile_1d = profile_3d(:,:,1);
        profile_all(:,offset) = Vprofile_1d;
        offset = offset + 1;
        for i=0:1:(prof_int_num-1)
            x0 = x0 + prof_int*catx/len;
            y0 = y0 - prof_int*caty/len;
            x1 = x1 + prof_int*catx/len;
            y1 = y1 - prof_int*caty/len;
            xi= [x0 x1];
            yi= [y0 y1];
            line(xi,yi);
            [cx, cy, profile_3d] = improfile(I,xi,yi,fp);
            Vprofile_1d = profile_3d(:,:,1);
            profile_all(:,offset) = Vprofile_1d;        
            offset = offset + 1;
        end   
    elseif gline_coordinates(1)>gline_coordinates(2) && gline_coordinates(3)>gline_coordinates(4)
        [pl_1, pl_2, pl_3, pl_4, xi, yi, x0, y0, x1, y1] = lines(gline_coordinates, catx, caty, hp, len); % See functions at the end of the file        
        confirm_msg = 'Is the profile line of the right length: ';
        confirm = strcat(confirm_msg, sprintf(' %d pixels (y/n)? ', hp));
        answer_1 = input(confirm, 's');
        while answer_1 == 'n'     
            new_length = 'Enter the half profile length? ';
            hp = input(new_length);
            fp = 2*hp;
            delete(pl_1);
            delete(pl_2);
            delete(pl_3);
            delete(pl_4);
            [pl_1, pl_2, pl_3, pl_4, xi, yi, x0, y0, x1, y1] = lines(gline_coordinates, catx, caty, hp, len);
            confirm_msg_1 = 'Is the profile line of the right length: ';
            confirm_1 = strcat(confirm_msg_1, sprintf(' %d pixels (y/n)? ', hp));
            answer_1 = input(confirm_1, 's'); 
        end
        delete(pl_2);
        delete(pl_3);
        delete(pl_4);        
        profile_all = zeros(fp, prof_int_num_temp); % Create empty matrix for all the profile values
        [cx, cy, profile_3d] = improfile(I,xi,yi,fp);
        hold on;
        plot(x0, y0, 'Marker','o','MarkerFaceColor','red');
        hold off;
        Vprofile_1d = profile_3d(:,:,1);
        profile_all(:,offset) = Vprofile_1d;
        offset = offset + 1;
        for i=0:1:(prof_int_num-1)
            x0 = x0 - prof_int*catx/len;
            y0 = y0 - prof_int*caty/len;
            x1 = x1 - prof_int*catx/len;
            y1 = y1 - prof_int*caty/len;
            xi= [x0 x1];
            yi= [y0 y1];
            line(xi,yi);
            [cx, cy, profile_3d] = improfile(I,xi,yi,fp);
            Vprofile_1d = profile_3d(:,:,1);
            profile_all(:,offset) = Vprofile_1d;        
            offset = offset + 1;
        end
    elseif gline_coordinates(1)<gline_coordinates(2) && gline_coordinates(3)<gline_coordinates(4)
        [pl_1, pl_2, pl_3, pl_4, xi, yi, x0, y0, x1, y1] = lines(gline_coordinates, catx, caty, hp, len); % See functions at the end of the file        
        confirm_msg = 'Is the profile line of the right length: ';
        confirm = strcat(confirm_msg, sprintf(' %d pixels (y/n)? ', hp));
        answer_1 = input(confirm, 's');
        while answer_1 == 'n'     
            new_length = 'Enter the half profile length? ';
            hp = input(new_length);
            fp = 2*hp;
            delete(pl_1);
            delete(pl_2);
            delete(pl_3);
            delete(pl_4);
            [pl_1, pl_2, pl_3, pl_4, xi, yi, x0, y0, x1, y1] = lines(gline_coordinates, catx, caty, hp, len);
            confirm_msg_1 = 'Is the profile line of the right length: ';
            confirm_1 = strcat(confirm_msg_1, sprintf(' %d pixels (y/n)? ', hp));
            answer_1 = input(confirm_1, 's'); 
        end
        delete(pl_2);
        delete(pl_3);
        delete(pl_4);        
        profile_all = zeros(fp, prof_int_num_temp); % Create empty matrix for all the profile values
        [cx, cy, profile_3d] = improfile(I,xi,yi,fp);
        hold on;
        plot(x0, y0, 'Marker','o','MarkerFaceColor','red');
        hold off;
        Vprofile_1d = profile_3d(:,:,1);
        profile_all(:,offset) = Vprofile_1d;
        offset = offset + 1;
        for i=0:1:(prof_int_num-1)
            x0 = x0 + prof_int*catx/len;
            y0 = y0 + prof_int*caty/len;
            x1 = x1 + prof_int*catx/len;
            y1 = y1 + prof_int*caty/len;
            xi= [x0 x1];
            yi= [y0 y1];
            line(xi,yi);
            [cx, cy, profile_3d] = improfile(I,xi,yi,fp);
            Vprofile_1d = profile_3d(:,:,1);
            profile_all(:,offset) = Vprofile_1d;        
            offset = offset + 1;
        end   
    elseif gline_coordinates(1)>gline_coordinates(2) && gline_coordinates(3)<gline_coordinates(4)
        [pl_1, pl_2, pl_3, pl_4, xi, yi, x0, y0, x1, y1] = lines(gline_coordinates, catx, caty, hp, len); % See functions at the end of the file        
        confirm_msg = 'Is the profile line of the right length: ';
        confirm = strcat(confirm_msg, sprintf(' %d pixels (y/n)? ', hp));
        answer_1 = input(confirm, 's');
        while answer_1 == 'n'     
            new_length = 'Enter the half profile length? ';
            hp = input(new_length);
            fp = 2*hp;
            delete(pl_1);
            delete(pl_2);
            delete(pl_3);
            delete(pl_4);
            [pl_1, pl_2, pl_3, pl_4, xi, yi, x0, y0, x1, y1] = lines(gline_coordinates, catx, caty, hp, len);
            confirm_msg_1 = 'Is the profile line of the right length: ';
            confirm_1 = strcat(confirm_msg_1, sprintf(' %d pixels (y/n)? ', hp));
            answer_1 = input(confirm_1, 's'); 
        end
        delete(pl_2);
        delete(pl_3);
        delete(pl_4);        
        profile_all = zeros(fp, prof_int_num_temp); % Create empty matrix for all the profile values
        [cx, cy, profile_3d] = improfile(I,xi,yi,fp);
        hold on;
        plot(x0, y0, 'Marker','o','MarkerFaceColor','red');
        hold off;
        Vprofile_1d = profile_3d(:,:,1);
        profile_all(:,offset) = Vprofile_1d;
        offset = offset + 1;
        for i=0:1:(prof_int_num-1)
            x0 = x0 - prof_int*catx/len;
            y0 = y0 + prof_int*caty/len;
            x1 = x1 - prof_int*catx/len;
            y1 = y1 + prof_int*caty/len;
            xi= [x0 x1];
            yi= [y0 y1];
            line(xi,yi);
            [cx, cy, profile_3d] = improfile(I,xi,yi,fp);
            Vprofile_1d = profile_3d(:,:,1);
            profile_all(:,offset) = Vprofile_1d;        
            offset = offset + 1;
        end         
    % Vertical lines
    elseif gline_coordinates(1) == gline_coordinates(2) && gline_coordinates(3)>gline_coordinates(4)
        x0 = gline_coordinates(1) + hp;
        y0 = gline_coordinates(3);
        x1 = gline_coordinates(1) - hp;
        y1 = gline_coordinates(3);
        xi= [x0 x1];
        yi= [y0 y1];
        prof_line = line(xi,yi);
        confirm_msg = 'Is the profile line of the right length: ';
        confirm = strcat(confirm_msg, sprintf('%d pixels (y/n)? ', hp));
        answer_1 = input(confirm, 's');
        while answer_1 == 'n'     
            new_length = 'Enter the half profile length? ';
            hp = input(new_length);
            fp = 2*hp;
            x0 = gline_coordinates(1) + hp;
            y0 = gline_coordinates(3);
            x1 = gline_coordinates(1) - hp;
            y1 = gline_coordinates(3);
            xi= [x0 x1];
            yi= [y0 y1];
            delete(prof_line);            
            prof_line = line(xi,yi);
            confirm_msg_1 = 'Is the profile line of the right length: ';
            confirm_1 = strcat(confirm_msg_1, sprintf('%d pixels (y/n)? ', hp));
            answer_1 = input(confirm_1, 's'); 
        end
        profile_all = zeros(fp, prof_int_num_temp); % Create empty matrix for all the profile values
        [cx, cy, profile_3d] = improfile(I,xi,yi,fp);
        Vprofile_1d = profile_3d(:,:,1);
        profile_all(:,offset) = Vprofile_1d;
        offset = offset + 1;
        for i=0:1:(prof_int_num-1)
            y0 = y0 - prof_int;
            y1 = y1 - prof_int;
            xi= [x0 x1];
            yi= [y0 y1];
            line(xi,yi);
            [cx, cy, profile_3d] = improfile(I,xi,yi,fp);
            Vprofile_1d = profile_3d(:,:,1);
            profile_all(:,offset) = Vprofile_1d;
            offset = offset + 1;
        end
    elseif gline_coordinates(1) == gline_coordinates(2) && gline_coordinates(3)<gline_coordinates(4)
        x0 = gline_coordinates(1) - hp;
        y0 = gline_coordinates(3);
        x1 = gline_coordinates(1) + hp;
        y1 = gline_coordinates(3);
        xi= [x0 x1];
        yi= [y0 y1];
        prof_line = line(xi,yi);
        confirm_msg = 'Is the profile line of the right length: ';
        confirm = strcat(confirm_msg, sprintf('%d pixels (y/n)? ', hp));
        answer_1 = input(confirm, 's');
        while answer_1 == 'n'     
            new_length = 'Enter the half profile length? ';
            hp = input(new_length);
            fp = 2*hp;
            x0 = gline_coordinates(1) - hp;
            y0 = gline_coordinates(3);
            x1 = gline_coordinates(1) + hp;
            y1 = gline_coordinates(3);
            xi= [x0 x1];
            yi= [y0 y1];
            delete(prof_line);            
            prof_line = line(xi,yi);
            confirm_msg_1 = 'Is the profile line of the right length: ';
            confirm_1 = strcat(confirm_msg_1, sprintf('%d pixels (y/n)? ', hp));
            answer_1 = input(confirm_1, 's'); 
        end
        profile_all = zeros(fp, prof_int_num_temp); % Create empty matrix for all the profile values
        [cx, cy, profile_3d] = improfile(I,xi,yi,fp);
        Vprofile_1d = profile_3d(:,:,1);
        profile_all(:,offset) = Vprofile_1d;
        offset = offset + 1;
        for i=0:1:(prof_int_num-1)
            y0 = y0 + prof_int;
            y1 = y1 + prof_int;
            xi= [x0 x1];
            yi= [y0 y1];
            line(xi,yi);
            [cx, cy, profile_3d] = improfile(I,xi,yi,fp);
            Vprofile_1d = profile_3d(:,:,1);
            profile_all(:,offset) = Vprofile_1d;
            offset = offset + 1;
        end        
    % Horizontal lines
    elseif gline_coordinates(3)== gline_coordinates(4) && gline_coordinates(1)<gline_coordinates(2)
        x0 = gline_coordinates(1);
        y0 = gline_coordinates(3) + hp;
        x1 = gline_coordinates(1);
        y1 = gline_coordinates(3) - hp;
        xi= [x0 x1];
        yi= [y0 y1];
        prof_line = line(xi,yi);
        confirm_msg = 'Is the profile line of the right length: ';
        confirm = strcat(confirm_msg, sprintf('%d pixels (y/n)? ', hp));
        answer_1 = input(confirm, 's');
        while answer_1 == 'n'     
            new_length = 'Enter the half profile length? ';
            hp = input(new_length);
            fp = 2*hp;
            x0 = gline_coordinates(1);
            y0 = gline_coordinates(3) + hp;
            x1 = gline_coordinates(1);
            y1 = gline_coordinates(3) - hp;
            xi= [x0 x1];
            yi= [y0 y1];
            delete(prof_line);            
            prof_line = line(xi,yi);
            confirm_msg_1 = 'Is the profile line of the right length: ';
            confirm_1 = strcat(confirm_msg_1, sprintf('%d pixels (y/n)? ', hp));
            answer_1 = input(confirm_1, 's'); 
        end
        profile_all = zeros(fp, prof_int_num_temp); % Create empty matrix for all the profile values
        [cx, cy, profile_3d] = improfile(I,xi,yi,fp);
        Vprofile_1d = profile_3d(:,:,1);
        profile_all(:,offset) = Vprofile_1d;
        offset = offset + 1;
        for i=0:1:(prof_int_num-1)
            x0 = x0 + prof_int;
            x1 = x1 + prof_int;
            xi= [x1 x0];
            yi= [y1 y0];
            line(xi,yi);
            [cx, cy, profile_3d] = improfile(I,xi,yi,fp);
            Vprofile_1d = profile_3d(:,:,1);
            profile_all(:,offset) = Vprofile_1d;
            offset = offset + 1;
        end
     elseif gline_coordinates(3)== gline_coordinates(4) && gline_coordinates(1)>gline_coordinates(2)
        x0 = gline_coordinates(1);
        y0 = gline_coordinates(3) - hp;
        x1 = gline_coordinates(1);
        xi= [x0 x1];
        yi= [y0 y1];
        prof_line = line(xi,yi);
        confirm_msg = 'Is the profile line of the right length: ';
        confirm = strcat(confirm_msg, sprintf('%d pixels (y/n)? ', hp));
        answer_1 = input(confirm, 's');
        while answer_1 == 'n'     
            new_length = 'Enter the half profile length? ';
            hp = input(new_length);
            fp = 2*hp;
            x0 = gline_coordinates(1);
            y0 = gline_coordinates(3) - hp;
            x1 = gline_coordinates(1);
            y1 = gline_coordinates(3) + hp;
            xi= [x0 x1];
            yi= [y0 y1];
            delete(prof_line);            
            prof_line = line(xi,yi);
            confirm_msg_1 = 'Is the profile line of the right length: ';
            confirm_1 = strcat(confirm_msg_1, sprintf('%d pixels (y/n)? ', hp));
            answer_1 = input(confirm_1, 's'); 
        end
        profile_all = zeros(fp, prof_int_num_temp); % Create empty matrix for all the profile values
        [cx, cy, profile_3d] = improfile(I,xi,yi,fp);
        Vprofile_1d = profile_3d(:,:,1);
        profile_all(:,offset) = Vprofile_1d;
        offset = offset + 1;
        for i=0:1:(prof_int_num-1)
            x0 = x0 - prof_int;
            x1 = x1 - prof_int;
            xi= [x0 x1];
            yi= [y0 y1];
            line(xi,yi);
            [cx, cy, profile_3d] = improfile(I,xi,yi,fp);
            Vprofile_1d = profile_3d(:,:,1);
            profile_all(:,offset) = Vprofile_1d;
            offset = offset + 1;
        end        
    end
    
    % Calculate the pixel distance between grey values and create an empty matrix to store real distances
    points_coordinates = horzcat(cx,cy);
    points_distance = pdist(points_coordinates);
    dist_pixel = zeros(fp, 1);
    for i=1:1:(fp-1)
        dist_pixel(i+1) = points_distance(i);
    end
    dist_real = zeros(fp, 1); % Create empty matrix for the real distances to be added
    profile_points = repmat(double(prof_int_num_temp),[fp,1]);
    
    % Perform calculations over grey values 
    mean_val = mean(profile_all,2);  
    min_val = min(profile_all,[],2);
    max_val = max(profile_all,[],2);
    std_dev = std(profile_all,0,2); 
    std_rel = std_dev ./ mean_val;
    std_error = std_dev ./ sqrt(profile_points);
    data_all_start = horzcat(profile_all, dist_pixel, dist_real, profile_points, profile_points, min_val, max_val, mean_val, std_dev, std_rel, std_error);
    
    % Exclude outliers and re-perform calculations 
    profile_all_calc = profile_all; % keep original values
% Uncomment the following lines to exclude grey values using digit comparison rather than difference.
%{    
    digits = 4;
    std_rel_str = cell(length(std_rel),1);        
    for i=1:1:length(std_rel_str)
        std_rel_str(i) = cellstr(num2str(std_rel(i)));
    end
    for i=1:1:(length(std_rel_str)-1)
        check_std = strncmp(std_rel_str(i),std_rel_str(i+1),digits);
        if check_std == 0 
            clean = true;
            break
        end
    end
%}  
    used_points = zeros(fp,1);
    nan_count = 0;
    for i=1:1:fp
        for j=1:1:(prof_int_num_temp)
            if profile_all_calc(i,j) > (mean_val(i) + std_dev(i)) || profile_all_calc(i,j) < (mean_val(i) - std_dev(i)) && profile_all_calc(i,j) ~= 0
                profile_all_calc(i,j) = NaN;
                nan_count = nan_count + 1;
            end
        end
        used_points(i) = prof_int_num_temp - nan_count;
        nan_count = 0;
    end
    mean_val = nanmean(profile_all_calc,2);  
    min_val = nanmin(profile_all_calc,[],2);
    max_val = nanmax(profile_all_calc,[],2);
    std_dev = nanstd(profile_all_calc,0,2); 
    std_rel = std_dev ./ mean_val;
    std_error = std_dev ./ sqrt(used_points);
    
    % Create a matrix with all the new values combined to export 
    data_all = horzcat(profile_all_calc, dist_pixel, dist_real, profile_points, used_points, min_val, max_val, mean_val, std_dev, std_rel, std_error);
      
    % Create a matrix with the titles for the CSV files 
    line_num = 1;
    line_array = cell(1, prof_int_num_temp);
    for l=1:1:(prof_int_num_temp)
        label = sprintf('Line_%d', line_num);
        line_array(:,l) = {label};
        line_num = line_num + 1;
    end
    label_rest = cell(1,10);
    label_rest(1,1) = {'Pixel_Distance'};
    label_rest(1,2) = {'Real_Distance'};
    label_rest(1,3) = {'Profiles'};
    label_rest(1,4) = {'Profiles_used_for_calculations'};
    label_rest(1,5) = {'Min_Value'};
    label_rest(1,6) = {'Max_Value'};
    label_rest(1,7) = {'Mean_Value'};
    label_rest(1,8) = {'Standard_Deviation'};
    label_rest(1,9) = {'Relative_Standard_Deviation'};
    label_rest(1,10) = {'Standard_Error'};
    label_all = horzcat(line_array, label_rest);     
    
    % Make a new directory checking if directory name exists, if it does, change new directory's name
    dir_number = 1; % Assign a number to the directory to avoid duplicates
    directory_name = strcat(ImageName(1:end-4), sprintf('_profile%d', dir_number));
    while (exist(directory_name, 'dir') == 7) 
        dir_number = dir_number + 1;
        directory_name = strcat(ImageName(1:end-4), sprintf('_profile%d', dir_number));
    end    
    name = strcat(ImageName(1:end-4), sprintf('_profile%d.csv', dir_number));
    name2 = strcat(ImageName(1:end-4), sprintf('_profile%d_clean.csv', dir_number));
    directory_name = name(1:end-4);
    mkdir(directory_name);
    directory_path = strcat('./', directory_name,'/');
    
    % Export titles and values to CSV files
    cell2csv(name,label_all); % See function at the end of the file. Credits to Sylvain Fiedler and Jerry Zhu
    cell2csv(name2,label_all);
    dlmwrite(name, data_all_start, '-append','delimiter',',');
    dlmwrite(name2, data_all, '-append','delimiter',',');
    
    % Create an image with all the profile lines
    name_fig = strcat(name(1:end-4),'.jpg');
    h = gcf;
    print (h, '-djpeg', name_fig);
    
    % Export guideline coordineates as matlab workspace
    full_gline_name = strcat(name(1:end-4),'_gline_coor.mat');
    save(full_gline_name, 'gline_coordinates');      
    
    % Move all files to the directory
    movefile(name, directory_path);
    movefile(name2, directory_path);
    movefile(full_gline_name, directory_path);
    movefile(name_fig, directory_path);   
    
    % Calculate time taken by the code to run and display it (useful for programming purposes)
    % time_taken = toc;
    % time_string = sprintf('Time taken for the code to run: %f Seconds', time_taken);
    % disp(time_string);
    
    % Ask what to do next    
    choice = menu('What Next?','Calculate Real World Distance (NB: scale on the image must be in microns)','Calculate more grey values', 'Exit The Program');
    
    % Calculate real world distances 
    % NB: this is based on the scale line usually present on the image used to
    % extrapolate greyvalues.
    while choice == 1
        [DistName] = uigetfile({'*.jpg;*.bmp;*.gif;*.png;*.tif','Images (*.jpg,*.bmp,*.gif,*.png,*.tif)'},'Select the IMAGE');
        imshow(DistName);
        figure(gcf); 
        % Uncomment to use ginput instead of imline
        % [x_dist_guideline,y_dist_guideline] = ginput(2);
        % line(x_dist_guideline,y_dist_guideline);
        % gline_dist_coordinates = horzcat(x_dist_guideline,y_dist_guideline);
        f = imline;
        gline_dist_coordinates = wait(f);
        distanceInPixels = sqrt( (gline_dist_coordinates(1)-gline_dist_coordinates(2)).^2 + (gline_dist_coordinates(3)-gline_dist_coordinates(4)).^2);

        % Ask the user for the real-world distance corresponding to the
        % scale line. Credits to Image Analyst. See function calibrate() at https://it.mathworks.com/matlabcentral/answers/uploaded_files/3815/spatial_calibration_demo.m 
        userPrompt = {'Enter real world units (e.g. microns):','Enter distance in those units:'};
        dialogTitle = 'Specify calibration information';
        numberOfLines = 1;
        def = {'microns', '500'};
        answer = inputdlg(userPrompt, dialogTitle, numberOfLines, def);
        if isempty(answer)
            return;
        end
        units = answer{1};
        distanceInUnits = str2double(answer{2});
        distancePerPixel = distanceInUnits / distanceInPixels;
        message = sprintf('The distance you drew is %.2f pixels = %f %s.\nThe number of %s per pixel is %f.\nThe number of pixels per %s is %f',...
        distanceInPixels, distanceInUnits, units, ...
        units, distancePerPixel, ...
        units, 1/distancePerPixel);
        uiwait(msgbox(message));
        
        % Write real distance data in the CSV file 
        dist_real = distancePerPixel*dist_pixel;
        data_all = horzcat(profile_all_calc, dist_pixel, dist_real, profile_points, used_points, min_val, max_val, mean_val, std_dev,std_rel,std_error);
        dist_real_column_number = prof_int_num_temp + 2;
        data_all_start(:,dist_real_column_number) = dist_real;
        cell2csv(name,label_all);
        cell2csv(name2,label_all);
        dlmwrite(name, data_all_start, '-append','delimiter',',');
        dlmwrite(name2, data_all, '-append','delimiter',',');
        
        % Delete old files and move new files to the directory
        filepath = strcat(directory_path, name);
        delete(filepath);
        filepath = strcat(directory_path, name2);
        delete(filepath);
        movefile(name, directory_path);
        movefile(name2, directory_path);
        choice = menu('What Next?','Calculate Real World Distance','Calculate more grey values', 'Exit The Program');
    end
    
    % Exit the program    
    if choice == 3 || choice == 0
        return;
    end 
end
end
   

%=================================================================================================================
%=================================================================================================================
%=================================================================================================================

% NB imread: y increaes downwards, x increases rightwards 

% Draws border lines do delineate the entire area over which grey profiles will be extracted.
function [pf_1, pf_2, pf_3, pf_4, Xi_1, Yi_1, X0, Y0, X1, Y1] = lines(position, cat1, cat2, half, length)

if position(1)<position(2) && position(3)>position(4) 
    X0 = position(1) + cat2*half/length; 
    Y0 = position(3) + cat1*half/length; 
    X1 = position(1) - cat2*half/length; 
    Y1 = position(3) - cat1*half/length; 
    X2 = position(2) + cat2*half/length; 
    Y2 = position(4) + cat1*half/length; 
    X3 = position(2) - cat2*half/length; 
    Y3 = position(4) - cat1*half/length; 
elseif  position(1)<position(2) && position(3)<position(4) 
    X0 = position(1) - cat2*half/length; 
    Y0 = position(3) + cat1*half/length; 
    X1 = position(1) + cat2*half/length; 
    Y1 = position(3) - cat1*half/length; 
    X2 = position(2) - cat2*half/length; 
    Y2 = position(4) + cat1*half/length; 
    X3 = position(2) + cat2*half/length; 
    Y3 = position(4) - cat1*half/length; 
elseif position(1)>position(2) && position(3)<position(4) 
    X0 = position(1) - cat2*half/length; 
    Y0 = position(3) - cat1*half/length; 
    X1 = position(1) + cat2*half/length; 
    Y1 = position(3) + cat1*half/length; 
    X2 = position(2) - cat2*half/length; 
    Y2 = position(4) - cat1*half/length; 
    X3 = position(2) + cat2*half/length; 
    Y3 = position(4) + cat1*half/length;  
elseif position(1)>position(2) && position(3)>position(4) 
    X0 = position(1) + cat2*half/length; 
    Y0 = position(3) - cat1*half/length; 
    X1 = position(1) - cat2*half/length; 
    Y1 = position(3) + cat1*half/length; 
    X2 = position(2) + cat2*half/length; 
    Y2 = position(4) - cat1*half/length; 
    X3 = position(2) - cat2*half/length; 
    Y3 = position(4) + cat1*half/length;    
end

Xi_1= [X0 X1];
Yi_1= [Y0 Y1];
Xi_2= [X2 X3];
Yi_2= [Y2 Y3];
Xi_3= [X0 X2];
Yi_3= [Y0 Y2];
Xi_4= [X1 X3];
Yi_4= [Y1 Y3];        
pf_1 = line(Xi_1,Yi_1, 'LineWidth',1.5,'Color', 'r');
pf_2 = line(Xi_2,Yi_2);
pf_3 = line(Xi_3,Yi_3);
pf_4 = line(Xi_4,Yi_4);

end

%=================================================================================================================
%=================================================================================================================
%=================================================================================================================

% Credits to Sylvain Fiedler and Jerry Zhu
% http://www.mathworks.com/matlabcentral/fileexchange/47055-cell-array-to-csv-file--improved-cell2csv-m-/content//cell2csv.m

function cell2csv(fileName, cellArray, separator, excelYear, decimal)

% % Writes cell array content into a *.csv file.
% % 
% % CELL2CSV(fileName, cellArray[, separator, excelYear, decimal])
% %
% % fileName     = Name of the file to save. [ e.g. 'text.csv' ]
% % cellArray    = Name of the Cell Array where the data is in
% % 
% % optional:
% % separator    = sign separating the values (default = ',')
% % excelYear    = depending on the Excel version, the cells are put into
% %                quotes before they are written to the file. The separator
% %                is set to semicolon (;)  (default = 1997 which does not change separator to semicolon ;)
% % decimal      = defines the decimal separator (default = '.')
% %
% %         by Sylvain Fiedler, KA, 2004
% % updated by Sylvain Fiedler, Metz, 06
% % fixed the logical-bug, Kaiserslautern, 06/2008, S.Fiedler
% % added the choice of decimal separator, 11/2010, S.Fiedler
% % modfiedy and optimized by Jerry Zhu, June, 2014, jerryzhujian9@gmail.com
% % now works with empty cells, numeric, char, string, row vector, and logical cells. 
% % row vector such as [1 2 3] will be separated by two spaces, that is "1  2  3"
% % One array can contain all of them, but only one value per cell.
% % 2x times faster than Sylvain's codes (8.8s vs. 17.2s):
% % tic;C={'te','tm';5,[1,2];true,{}};C=repmat(C,[10000,1]);cell2csv([datestr(now,'MMSS') '.csv'],C);toc;

%% Checking for optional Variables
if ~exist('separator', 'var')
    separator = ',';
end

if ~exist('excelYear', 'var')
    excelYear = 1997;
end

if ~exist('decimal', 'var')
    decimal = '.';
end

%% Setting separator for newer excelYears
if excelYear > 2000
    separator = ';';
end

% convert cell
cellArray = cellfun(@StringX, cellArray, 'UniformOutput', false);

%% Write file
datei = fopen(fileName, 'w');
[nrows,ncols] = size(cellArray);
for row = 1:nrows
    fprintf(datei,[sprintf(['%s' separator],cellArray{row,1:ncols-1}) cellArray{row,ncols} '\n']);
end    
% Closing file
fclose(datei);

% sub-function
function x = StringX(x)
    % If zero, then empty cell
    if isempty(x)
        x = '';
    % If numeric -> String, e.g. 1, [1 2]
    elseif isnumeric(x) && isrow(x)
        x = num2str(x);
        if decimal ~= '.'
            x = strrep(x, '.', decimal);
        end
    % If logical -> 'true' or 'false'
    elseif islogical(x)
        if x == 1
            x = 'TRUE';
        else
            x = 'FALSE';
        end
    % If matrix array -> a1 a2 a3. e.g. [1 2 3]
    % also catch string or char here
    elseif isrow(x) && ~iscell(x)
        x = num2str(x);
    % everthing else, such as [1;2], {1}
    else
        x = 'NA';
    end

    % If newer version of Excel -> Quotes 4 Strings
    if excelYear > 2000
        x = ['"' x '"'];
    end
end % end sub-function
end % end function


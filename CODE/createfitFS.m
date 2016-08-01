% This code uses MATLAB Curve Fitting TOOLBOX to fit a slightly modified
% error function to sample points. Created to be used with greyvalues
% and calculate the residence time of crystals in magmatic chambers.

% Residence time is calculated using given temperature, diffusion coefficient and
% pressure. These values depend on the mineral. Default values are for pyroxene 
% crystals. The error on residence time is calculated using different
% methods (e.g.: error propagation, normal and uniform distribution).

% N.B.: to be used with finite sources profiles (i.e.: gaussian shape
% profiles). It hence REQUIRES the input of a previously found parameter. 
% That is 'y0' found by sister function 'createfit'. 
% This parameter is used to calculate a value for 'erf_par' which is set fix in fitting the sample points. 

% Similarly, the code requires selecting a CSV input file. This MUST be formatted in
% the same way as the output file from sister function 'greyvalues'.

%  Output folder:
%      Data on fitting parameters: y0, h, erf_par, x0, Dt. 
%      Data on residence time calculations.
%      Image of the fit.

% The code has been run on machines running Matlab 2011b, 2014b and 2015a and 
% Windows XP, Windows 7, Windows 10, Ubuntu 12.04 and OS 10.7 and 10.10. 


function [fitresult, gof] = createfitFS()

% Select CSV file with data to fit
[csvname,csvpath] = uigetfile('*.csv','Select the CSV file with greyvalues to fit');
data_to_fit = csvimport(strcat(csvpath,csvname)); % See function at the end of the file. Credits to Ashish Sadanandan

% Set up values from CSV file for the fit: distance, greyvalues and relative errors
for i=1:1:size(data_to_fit,2)
    if strcmp(data_to_fit{1,i},'Real_Distance') == 1
        distance_cell = data_to_fit(:,i);
        distance = zeros((size(data_to_fit,1)-1),1);
        for k=2:1:size(data_to_fit,1)
            distance(k-1,1) = cell2mat(distance_cell(k,1));
        end
    elseif strcmp(data_to_fit{1,i},'Mean_Value') == 1
        greyvalues_cell = data_to_fit(:,i);
        greyvalues = zeros((size(data_to_fit,1)-1),1);
        for k=2:1:size(data_to_fit,1)
            greyvalues(k-1,1) = cell2mat(greyvalues_cell(k,1));        
        end
    elseif strcmp(data_to_fit{1,i},'Standard_Error') == 1
        errors_cell = data_to_fit(:,i);
        errors = zeros((size(data_to_fit,1)-1),1);
        for k=2:1:size(data_to_fit,1)
            errors(k-1,1) = cell2mat(errors_cell(k,1));
        end        
    end
end

% Set up starting data and initial guesses for parameters 
weights = 1 ./ (errors.^2); 
xdata = distance;
ydata = greyvalues;
middle_index = round(length(xdata)/2);
x0_start = xdata(middle_index);
y0_start = ydata(1);  
h_start = 5;
erf_param_start_msg = 'Enter initial condition (y0 previously modelled): ';
y0_old = input(erf_param_start_msg);
erf_par = abs(y0_old - ydata(1)) / 2; % Fixed value based on WRONG 'y0_start' estimate.
Dt_start = 6;

% Set up fittype and options to get the fit value of y0, used to recalculate 'erf_par' 
equation_ERFC = sprintf('y0+(%d*(erfc((h-(x-x0))/Dt) + erfc((h+(x-x0))/Dt)))',erf_par);
equation_ERF = sprintf('y0+(%d*(erf((h-(x-x0))/Dt) + erf((h+(x-x0))/Dt)))',erf_par);
if ydata(middle_index) < ydata(1) % ERFC 
    ft = fittype(equation_ERFC,...
    'dependent',{'y'},'independent',{'x'},...
    'coefficients',{'y0', 'h', 'x0', 'Dt'});
elseif ydata(middle_index) > ydata(1) %ERF
     ft = fittype(equation_ERF,...
    'dependent',{'y'},'independent',{'x'},...
    'coefficients',{'y0', 'h', 'x0', 'Dt'});
end
opts = fitoptions( ft );
opts.Weights = weights;
opts.DiffMaxChange = 0.01;
opts.Display = 'Off';
opts.MaxFunEvals = 1000000000;
opts.MaxIter = 1000000000;
opts.StartPoint = [y0_start h_start x0_start Dt_start];
opts.Upper = [Inf Inf Inf Inf];
opts.Lower = [-Inf -Inf -Inf -Inf];
% Fit model to data.
[fitresult, gof] = fit( xdata, ydata, ft, opts );
param_values = coeffvalues(fitresult);
erf_par = abs(y0_old - param_values(1)) / 2; % Fixed value based on RIGHT 'y0' value found by fit function 

% Re compute fitting using correct 'erf_par' 
equation_ERFC = sprintf('y0+(%d*(erfc((h-(x-x0))/Dt) + erfc((h+(x-x0))/Dt)))',erf_par);
equation_ERF = sprintf('y0+(%d*(erf((h-(x-x0))/Dt) + erf((h+(x-x0))/Dt)))',erf_par);
if ydata(middle_index) < ydata(1) % ERFC 
    ft = fittype(equation_ERFC,...
    'dependent',{'y'},'independent',{'x'},...
    'coefficients',{'y0', 'h', 'x0', 'Dt'});
elseif ydata(middle_index) > ydata(1) %ERF
     ft = fittype(equation_ERF,...
    'dependent',{'y'},'independent',{'x'},...
    'coefficients',{'y0', 'h', 'x0', 'Dt'});
end
% Fit model to data.
[fitresult, gof] = fit( xdata, ydata, ft, opts );

% Set up figure to receive data sets and fits
graph = clf;
img = figure(graph);
set(graph,'Units','Pixels','Position',[549 276 688 485]);
% Line handles and text for the legend.
legh = [];
legt = {};
% Limits of the x-axis.
xlim = [Inf -Inf];
% Axes for the plot.
ax = axes;
set(ax,'Units','normalized','OuterPosition',[0 0 1 1]);
set(ax,'Box','on');
axes(ax);
hold on;

% Plot data "Grey Values vs. Distance(um)"
xdata = xdata(:);
ydata = ydata(:);
graph_data = line(xdata,ydata,'Parent',ax,'Color',[0 0 1],...
    'linestyle', 'none',...
    'Marker','.', 'MarkerSize',15);
errorbar(xdata,ydata,errors, 'linestyle','none','Color',[0.5 0.5 0.5]);
xlim(1) = min(xlim(1),min(xdata));
xlim(2) = max(xlim(2),max(xdata));
legh(end+1) = graph_data;
legt{end+1} = 'Measured data';

% Nudge axis limits beyond data limits
if all(isfinite(xlim))
    xlim = xlim + [-1 1] * 0.01 * diff(xlim);
    set(ax,'XLim',xlim)
else
    set(ax, 'XLim',[-1, 101]);
end

% Plot the fit
fit_graph = plot(fitresult,'fit',0.95);
set(fit_graph(1),'Color',[1 0 0],...
    'LineStyle','-', 'LineWidth',2,...
    'Marker','none', 'MarkerSize',6);
legend off; % Turn off legend created by plot method.
% Store line handle and fit name for legend.
legh(end+1) = fit_graph(1);
legt{end+1} = 'Fit';

% Finish plotting data. 
hold off;
% Display legend
leginfo = {'Orientation', 'vertical', 'Location', 'NorthEast'};
fit_graph = legend(ax,legh,legt,leginfo{:});
set(fit_graph,'Interpreter','none');
% Label x- and y-axes.
xlabel(ax,'\mum');
ylabel(ax,'Greyvalues [0-255]');

% Make a new directory to export data on fitting parameters and residence time
dir_number = 1;
directory_name = strcat(csvname(1:end-4), sprintf('_fitFS_%d', dir_number));
while (exist(directory_name, 'dir') == 7)
    dir_number = dir_number + 1;
    directory_name = strcat(csvname(1:end-4), sprintf('_fitFS_%d', dir_number));
end
mkdir(directory_name);

% Export graph
img_name = sprintf('fit_graph_%d', dir_number);
img_name_full = strcat('./',directory_name, '/', img_name);
print(img, img_name_full, '-dpng');

% Export fitting data as variables to a workspace in the new directory.
conf_int = confint(fitresult);
param_names = coeffnames(fitresult);
param_values = coeffvalues(fitresult);
if ydata(middle_index) < ydata(1) % ERFC 
    yfit = param_values(1)+(erf_par*(erfc((param_values(2)-(xdata-param_values(3)))/param_values(4)) + erfc((param_values(2)+(xdata-param_values(3)))/param_values(4))));
elseif ydata(middle_index) > ydata(1) % ERF 
    yfit = param_values(1)+(erf_par*(erf((param_values(2)-(xdata-param_values(3)))/param_values(4)) + erf((param_values(2)+(xdata-param_values(3)))/param_values(4))));
end
y0_data = horzcat(param_names(1),param_values(1),conf_int(:,1));
erf_par_name = cellstr('ERF parameter set using input of initial conditions (i.e.: y0)');
erf_par_data = horzcat(erf_par_name,erf_par,sprintf('(%d - %d)/2 ',y0_old,param_values(1)));
h_data = horzcat(param_names(2),param_values(2),conf_int(:,2));
x0_data = horzcat(param_names(3),param_values(3),conf_int(:,3));
Dt_data = horzcat(param_names(4),param_values(4),conf_int(:,4));

workspace_name = sprintf('data_on_the_fit_%d.mat', dir_number);
workspace_full_name = strcat('./',directory_name, '/', workspace_name);
save(workspace_full_name, 'param_values', 'yfit', 'gof', 'x0_data', 'y0_data', 'h_data', 'erf_par_data', 'Dt_data');

% Export matlab variable with Dt value and confidence interval.
Dt_name = sprintf('Dt_data_%d.mat', dir_number);
Dt_full_name = strcat('./',directory_name, '/', Dt_name);
save(Dt_full_name, 'Dt_data');

% Calculate residence time and error on residence time using different methods
calc_msg = 'Would you like to calculate the residence time (D0 - Diffusion coefficient; DH - Pressure; Temperature and its range (DeltaT) are required?(y/n) ';
answer = input(calc_msg, 's');
if answer == 'n'
    return;
elseif answer == 'y'
    D0_msg = 'Use DO = 9.55E-05 m2s-1 ? (y/n) '; 
    D0_answer = input(D0_msg, 's');
    D0 = 9.55E-05;
    if D0_answer == 'n'
        D0_msg = 'Enter D0: ';
        D0 = input(D0_msg);
    end
    DH_msg = 'Use DeltaH = 406 kJmol-1 ? (y/n) '; 
    DH_answer = input(DH_msg, 's');
    DH = 406;
    if DH_answer == 'n'
        DH_msg = 'Enter DeltaH: ';
        DH = input(DH_msg);
    end
    T_msg = 'Enter Temperature (°C): ';
    t_c = input(T_msg);
    t_k = t_c + 273.15;
    
    DT_msg = 'Use DeltaT = 10°C ? (y/n) '; 
    DT_answer = input(DT_msg, 's');
    DeltaT = 10;
    if DT_answer == 'n'
        DT_msg = 'Enter DeltaT: ';
        DeltaT = input(DT_msg);
    end    
    
    % Dt was calculated using microns, turn it into a value in metres
    Dt = cell2mat(Dt_data(2))*0.000001;
    Dt_err = cell2mat(Dt_data(3))*0.000001;
    R =  8.31451/1000;
    
    % Maximum error
    D = D0*exp(-DH./(R*t_k));
    D_max = D0*exp(-DH./(R*(t_k-DeltaT)));
    D_err = (D - D_max) / D;
    res_time_s = (Dt/2)^2*(1./D);
    res_time_s_max = (Dt_err(2)/2)^2*(1/D_max);
    res_time_y = res_time_s/(60*60*24*365.25);
    res_time_max_y = res_time_s_max/(60*60*24*365.25);
    difference_max = (res_time_max_y - res_time_y)/res_time_y;

    % Use error propagation
    DDt = 2*(Dt_err(2)-Dt)/Dt;
    DD = DH/R * ((DeltaT/t_k)*1/t_k);
    t_err = (DDt^2+ DD^2)^0.5; % In percentage
    
    % Errors on Dt and temperature are uniformly distributed
    ran_val = 1000000;
    numbers_dt = (Dt_err(2)-Dt_err(1)).*rand(ran_val,1) + Dt_err(1);
    numbers = ((t_k+DeltaT)-(t_k-DeltaT)).*rand(ran_val,1) + (t_k-DeltaT);
    D_val = D0*exp(-DH./(R*numbers));
    D_val_err = std(D_val)/mean(D_val);
    res_time_s_val = numbers_dt.^2./(4*D_val);
    res_time_y_val = res_time_s_val./(60*60*24*365.25);
    res_time_val = mean(res_time_y_val);
    diff_val = std(res_time_y_val)/mean(res_time_y_val);
    
    % Errors on Dt and temperature are normally distributed
    nor_dt = Dt + ((Dt_err(2)-Dt).*randn(ran_val,1));
    nor_t = t_k + 10.*randn(ran_val,1);
    D_nor = D0*exp(-DH./(R*nor_t));
    D_nor_err = std(D_nor)/mean(D_nor);
    res_time_s_nor = nor_dt.^2./(4*D_nor);
    res_time_y_nor = res_time_s_nor./(60*60*24*365.25);
    res_time_nor = mean(res_time_y_nor);
    diff_nor = std(res_time_y_nor)/mean(res_time_y_nor);
    
    % Error on Dt is normally distributed. Error on temperature is uniformly distributed
    D_mix = D0*exp(-DH./(R*numbers));
    D_mix_err = std(D_mix)/mean(D_mix);
    res_time_s_mix = nor_dt.^2./(4*D_mix);
    res_time_y_mix = res_time_s_mix./(60*60*24*365.25);
    res_time_mix = mean(res_time_y_mix);
    diff_mix = std(res_time_y_mix)/mean(res_time_y_mix);
    
    % The temperature is constant. Error on Dt is normally distributed
    D_const = D0*exp(-DH./(R*t_k));
    D_const_err = 0;
    res_time_s_const = nor_dt.^2./(4*D_const);
    res_time_y_const = res_time_s_const./(60*60*24*365.25);
    res_time_const = mean(res_time_y_const);
    diff_const = std(res_time_y_const)/mean(res_time_y_const);
    
    % Export data on residence time and errors as variables to a workspace in the dedicated folder  
    blank = cellstr('----');
    temp_text = cellstr('Temperature_(°C)');
    res_time_text = cellstr('Residence_time_(y)');
    abs_error_text = cellstr('relative_error_on_residence_time(%)');
    rel_error_text = cellstr('relative_error_on_Dt_parameter(%)');
    key = vertcat(blank,temp_text,res_time_text,abs_error_text,rel_error_text);
    res_time_max_error_s = cellstr('res_time_max_error');
    res_time_max_error = vertcat(res_time_max_error_s,t_c,res_time_y,difference_max, D_err);
    res_time_error_prop_s = cellstr('res_time_error_prop');
    res_time_error_prop = vertcat(res_time_error_prop_s,t_c,res_time_y,t_err,DD);
    res_time_uniform_dist_s = cellstr('res_time_uniform_dist');
    res_time_uniform_dist = vertcat(res_time_uniform_dist_s,t_c,res_time_val,diff_val, D_val_err);
    res_time_normal_dist_s = cellstr('res_time_normal_dist');
    res_time_normal_dist = vertcat(res_time_normal_dist_s,t_c,res_time_nor,diff_nor, D_nor_err);
    res_time_temp_uniform_Dt_normal_s = cellstr('res_time_temp_uniform_Dt_normal');
    res_time_temp_uniform_Dt_normal = vertcat(res_time_temp_uniform_Dt_normal_s,t_c,res_time_mix,diff_mix, D_mix_err);
    res_time_temp_const_s = cellstr('res_time_temp_const');
    res_time_temp_const = vertcat(res_time_temp_const_s,t_c,res_time_y,DDt, D_const_err);
    res_time_temp_const_Dt_normal_s = cellstr('res_time_temp_const_Dt_normal');
    res_time_temp_const_Dt_normal = vertcat(res_time_temp_const_Dt_normal_s,t_c,res_time_const,diff_const, D_const_err);
    res_time_all = horzcat(key,res_time_max_error,res_time_error_prop,res_time_uniform_dist,res_time_normal_dist,res_time_temp_uniform_Dt_normal,res_time_temp_const,res_time_temp_const_Dt_normal);
    res_time_all_name = sprintf('res_time_all_%d.mat', dir_number);
    res_time_all_full_name = strcat('./',directory_name, '/', res_time_all_name);
    
    D0_s = cellstr('D0');
    D0_data = vertcat(D0_s,D0);
    DH_s = cellstr('DeltaH');
    DH_data = vertcat(DH_s,DH);
    T_s = cellstr('Temperature');
    T_data = vertcat(T_s,t_c);
    DT_s = cellstr('DeltaT');
    DT_data = vertcat(DT_s,DeltaT);
    data_input_all = horzcat(D0_data,DH_data,T_data,DT_data);

    save(res_time_all_full_name, 'res_time_all', 'data_input_all');
end    

% Compute calculations using different Dt and calculate error using only
% error propagation method.
yfit_inserted_Dt_msg = 'Would you like to calculate the fitting using a different Dt parameter (MUST select a .mat file with Dt data) ? (y/n) ';
yfit_inserted_Dt_msg_answer = input(yfit_inserted_Dt_msg, 's');
if yfit_inserted_Dt_msg_answer == 'n'
    return;
elseif yfit_inserted_Dt_msg_answer == 'y'
    [Dt_name,Dt_path] = uigetfile('*.mat','Select the file with Dt data to use');
    load(strcat(Dt_path,Dt_name)); 
    T_s = sprintf('Use temperature = %d °C ? (y/n)', t_c);
    if input(T_s, 's') == 'n'
        T_msg = 'Enter Temperature (in Celsius): ';
        t_c = input(T_msg);
        t_k = t_c + 273.15;
    end    
    % Dt was calculated using microns, turn it into a value in metres
    Dt = cell2mat(Dt_data(2))*0.000001;
    Dt_err = cell2mat(Dt_data(3))*0.000001;
    R =  8.31451/1000;
    D = D0*exp(-DH./(R*t_k));
    res_time_s = (Dt/2)^2*(1./D);
    res_time_y = res_time_s/(60*60*24*365.25);
    
    % Use error propagation
    DDt = 2*(Dt_err(2)-Dt)/Dt;
    DD = DH/R * ((DeltaT/t_k)*1/t_k);
    t_err = (DDt^2+ DD^2)^0.5; % In percentage
   
    % Export time data
    time_error_prop_s = cellstr('time_error_prop');
    time_error_prop = vertcat(time_error_prop_s,t_c,res_time_y,t_err, DD);
    time_tem_const_s = cellstr('time_tem_const');
    time_tem_const = vertcat(time_tem_const_s,t_c,res_time_y,DDt, D_const_err);
    data_time_all = horzcat(key,time_error_prop,time_tem_const);
   
    % Export new Dt parameters 
    Dt_inserted_Dt = Dt;
    param_values_inserted_Dt = param_values;
    param_values_inserted_Dt(4) = Dt_inserted_Dt;
    if ydata(middle_index) < ydata(1) % ERFC 
        yfit_inserted_Dt = param_values_inserted_Dt(1)+(erf_par*(erfc((param_values_inserted_Dt(2)-(xdata-param_values_inserted_Dt(3)))/param_values_inserted_Dt(4)) + erfc((param_values_inserted_Dt(2)+(xdata-param_values_inserted_Dt(3)))/param_values_inserted_Dt(4))));
    elseif ydata(middle_index) > ydata(1) % ERF 
        yfit_inserted_Dt = param_values_inserted_Dt(1)+(erf_par*(erf((param_values_inserted_Dt(2)-(xdata-param_values_inserted_Dt(3)))/param_values_inserted_Dt(4)) + erf((param_values_inserted_Dt(2)+(xdata-param_values_inserted_Dt(3)))/param_values_inserted_Dt(4))));
    end
    
    res_time_inserted_Dt_name = sprintf('res_time_inserted_Dt_%d.mat', dir_number);
    res_time_inserted_Dt_full_name = strcat('./',directory_name, '/', res_time_inserted_Dt_name);
    save(res_time_inserted_Dt_full_name, 'data_time_all', 'yfit_inserted_Dt', 'Dt_data', 'param_values_inserted_Dt'); 
end
end

%=================================================================================================================
%=================================================================================================================
%=================================================================================================================

% Credits to Ashish Sadanandan
% https://www.mathworks.com/matlabcentral/fileexchange/23573-csvimport

function varargout = csvimport( fileName, varargin )
% CSVIMPORT reads the specified CSV file and stores the contents in a cell array or matrix
%
% The file can contain any combination of text & numeric values. Output data format will vary
% depending on the exact composition of the file data.
%
% CSVIMPORT( fileName ):         fileName     -  String specifying the CSV file to be read. Set to
%                                                [] to interactively select the file.
%
% CSVIMPORT( fileName, ... ) : Specify a list of options to be applied when importing the CSV file.
%                              The possible options are:
%                                delimiter     - String to be used as column delimiter. Default
%                                                value is , (comma)
%                                columns       - String or cell array of strings listing the columns
%                                                from which data is to be extracted. If omitted data
%                                                from all columns in the file is imported. If file
%                                                does not contain a header row, the columns
%                                                parameter can be a numeric array listing column
%                                                indices from which data is to be extracted.
%                                outputAsChar  - true / false value indicating whether the data
%                                                should be output as characters. If set to false the
%                                                function attempts to convert each column into a
%                                                numeric array, it outputs the column as characters
%                                                if conversion of any data element in the column
%                                                fails. Default value is false.
%                                uniformOutput - true / false value indicating whether output can be
%                                                returned without encapsulation in a cell array.
%                                                This parameter is ignored if the columns / table
%                                                cannot be converted into a matrix.
%                                noHeader      - true / false value indicating whether the CSV
%                                                file's first line contains column headings. Default
%                                                value is false.
%                                ignoreWSpace  - true / false value indicating whether to ignore
%                                                leading and trailing whitespace in the column
%                                                headers; ignored if noHeader is set to true.
%                                                Default value is false.
%
% The parameters must be specified in the form of param-value pairs, parameter names are not
% case-sensitive and partial matching is supported.
%
% [C1 C2 C3] = CSVIMPORT( fileName, 'columns', {'C1', 'C2', C3'}, ... )
%   This form returns the data from columns in output variables C1, C2 and C3 respectively, the
%   column names are case-sensitive and must match a column name in the file exactly. When fetching
%   data in column mode the number of output columns must match the number of columns to read or it
%   must be one. In the latter case the data from the columns is returned as a single cell matrix.
%
% [C1 C2 C3] = CSVIMPORT( fileName, 'columns', [1, 3, 4], ,'noHeader', true, ... )
%   This form returns the data from columns in output variables C1, C2 and C3 respectively, the
%   columns parameter must contain the column indices when the 'noHeader' option is set to true.

%
% Notes:  1. Function has not been tested on badly formatted CSV files.
%         2. Created using R2007b but has been tested on R2006b.
%
% Revisions:
%   04/28/2009: Corrected typo in an error message
%               Added igonoreWSpace option
%   08/16/2010: Replaced calls to str2num with str2double, the former uses eval leading to unwanted
%               side effects if cells contain text with function names
%

if ( nargin == 0 ) || isempty( fileName )
  [fileName filePath] = uigetfile( '*.csv', 'Select CSV file' );
  if isequal( fileName, 0 )
    return;
  end
  fileName = fullfile( filePath, fileName );
else
  if ~ischar( fileName )
    error( 'csvimport:FileNameError', 'The first argument to %s must be a valid .csv file', ...
      mfilename );
  end
end

%Setup default values
p.delimiter       = ',';
p.columns         = [];
p.outputAsChar    = false;
p.uniformOutput   = true;
p.noHeader        = false;
p.ignoreWSpace    = false;

validParams     = {     ...
  'delimiter',          ...
  'columns',            ...
  'outputAsChar',       ...
  'uniformOutput',      ...
  'noHeader',           ...
  'ignoreWSpace'        ...
  };

%Parse input arguments
if nargin > 1
  if mod( numel( varargin ), 2 ) ~= 0
    error( 'csvimport:InvalidInput', ['All input parameters after the fileName must be in the ' ...
      'form of param-value pairs'] );
  end
  params  = lower( varargin(1:2:end) );
  values  = varargin(2:2:end);

  if ~all( cellfun( @ischar, params ) )
    error( 'csvimport:InvalidInput', ['All input parameters after the fileName must be in the ' ...
      'form of param-value pairs'] );
  end

  lcValidParams   = lower( validParams );
  for ii =  1 : numel( params )
    result        = strmatch( params{ii}, lcValidParams );
    %If unknown param is entered ignore it
    if isempty( result )
      continue
    end
    %If we have multiple matches make sure we don't have a single unambiguous match before throwing
    %an error
    if numel( result ) > 1
      exresult    = strmatch( params{ii}, validParams, 'exact' );
      if ~isempty( exresult )
        result    = exresult;
      else
        %We have multiple possible matches, prompt user to provide an unambiguous match
        error( 'csvimport:InvalidInput', 'Cannot find unambiguous match for parameter ''%s''', ...
          varargin{ii*2-1} );
      end
    end
    result      = validParams{result};
    p.(result)  = values{ii};
  end
end

%Check value attributes
if isempty( p.delimiter ) || ~ischar( p.delimiter )
  error( 'csvimport:InvalidParamType', ['The ''delimiter'' parameter must be a non-empty ' ...
    'character array'] );
end
if isempty( p.noHeader ) || ~islogical( p.noHeader ) || ~isscalar( p.noHeader )
  error( 'csvimport:InvalidParamType', ['The ''noHeader'' parameter must be a non-empty ' ...
    'logical scalar'] );
end
if ~p.noHeader
  if ~isempty( p.columns )
    if ~ischar( p.columns ) && ~iscellstr( p.columns )
      error( 'csvimport:InvalidParamType', ['The ''columns'' parameter must be a character array ' ...
        'or a cell array of strings for CSV files containing column headers on the first line'] );
    end
    if p.ignoreWSpace
      p.columns = strtrim( p.columns );
    end
  end
else
  if ~isempty( p.columns ) && ~isnumeric( p.columns )
    error( 'csvimport:InvalidParamType', ['The ''columns'' parameter must be a numeric array ' ...
      'for CSV files containing column headers on the first line'] );
  end
end
if isempty( p.outputAsChar ) || ~islogical( p.outputAsChar ) || ~isscalar( p.outputAsChar )
  error( 'csvimport:InvalidParamType', ['The ''outputAsChar'' parameter must be a non-empty ' ...
    'logical scalar'] );
end
if isempty( p.uniformOutput ) || ~islogical( p.uniformOutput ) || ~isscalar( p.uniformOutput )
  error( 'csvimport:InvalidParamType', ['The ''uniformOutput'' parameter must be a non-empty ' ...
    'logical scalar'] );
end

%Open file
[fid msg] = fopen( fileName, 'rt' );
if fid == -1
  error( 'csvimport:FileReadError', 'Failed to open ''%s'' for reading.\nError Message: %s', ...
    fileName, msg );
end

colMode         = ~isempty( p.columns );
if ischar( p.columns )
  p.columns     = cellstr( p.columns );
end
nHeaders        = numel( p.columns );

if colMode
  if ( nargout > 1 ) && ( nargout ~= nHeaders )
    error( 'csvimport:NumOutputs', ['The number of output arguments must be 1 or equal to the ' ...
      'number of column names when fetching data for specific columns'] );
  end
end

%Read first line and determine number of columns in data
rowData         = fgetl( fid );
rowData         = regexp( rowData, p.delimiter, 'split' );
nCols           = numel( rowData );

%Check whether all specified columns are present if used in column mode and store their indices
if colMode
  if ~p.noHeader
    if p.ignoreWSpace
      rowData     = strtrim( rowData );
    end
    colIdx        = zeros( 1, nHeaders );
    for ii = 1 : nHeaders
      result      = strmatch( p.columns{ii}, rowData );
      if isempty( result )
        fclose( fid );
        error( 'csvimport:UnknownHeader', ['Cannot locate column header ''%s'' in the file ' ...
          '''%s''. Column header names are case sensitive.'], p.columns{ii}, fileName );
      elseif numel( result ) > 1
        exresult  = strmatch( p.columns{ii}, rowData, 'exact' );
        if numel( exresult ) == 1
          result  = exresult;
        else
          warning( 'csvimport:MultipleHeaderMatches', ['Column header name ''%s'' matched ' ...
            'multiple columns in the file, only the first match (C:%d) will be used.'], ...
            p.columns{ii}, result(1) );
        end
      end
      colIdx(ii)  = result(1);
    end
  else
    colIdx        = p.columns(:);
    if max( colIdx ) > nCols
      fclose( fid );
      error( 'csvimport:BadIndex', ['The specified column index ''%d'' exceeds the number of ' ...
        'columns (%d) in the file'], max( colIdx ), nCols );
    end
  end
end

%Calculate number of lines
pos             = ftell( fid );
if pos == -1
  msg = ferror( fid );
  fclose( fid );
  error( 'csvimport:FileQueryError', 'FTELL on file ''%s'' failed.\nError Message: %s', ...
    fileName, msg );
end
data            = fread( fid );
nLines          = numel( find( data == sprintf( '\n' ) ) ) + 1;
%Reposition file position indicator to beginning of second line
if fseek( fid, pos, 'bof' ) ~= 0
  msg = ferror( fid );
  fclose( fid );
  error( 'csvimport:FileSeekError', 'FSEEK on file ''%s'' failed.\nError Message: %s', ...
    fileName, msg );
end

data            = cell( nLines, nCols );
data(1,:)       = rowData;
emptyRowsIdx    = [];
%Get data for remaining rows
for ii = 2 : nLines
  rowData       = fgetl( fid );
  if isempty( rowData )
    emptyRowsIdx = [emptyRowsIdx(:); ii];
    continue
  end
  rowData       = regexp( rowData, p.delimiter, 'split' );
  nDataElems    = numel( rowData );
  if nDataElems < nCols
    warning( 'csvimport:UnevenColumns', ['Number of data elements on line %d (%d) differs from ' ...
      'that on the first line (%d). Data in this line will be padded.'], ii, nDataElems, nCols );
    rowData(nDataElems+1:nCols) = {''};
  elseif nDataElems > nCols
    warning( 'csvimport:UnevenColumns', ['Number of data elements on line %d (%d) differs from ' ...
      'that one the first line (%d). Data in this line will be truncated.'], ii, nDataElems, nCols );
    rowData     = rowData(1:nCols);
  end
  data(ii,:)    = rowData;
end
%Close file handle
fclose( fid );
data(emptyRowsIdx,:)   = [];

%Process data for final output
uniformOutputPossible  = ~p.outputAsChar;
if p.noHeader
  startRowIdx          = 1;
else
  startRowIdx          = 2;
end
if ~colMode
  if ~p.outputAsChar
    %If we're not outputting the data as characters then try to convert each column to a number
    for ii = 1 : nCols
      colData     = cellfun( @str2double, data(startRowIdx:end,ii), 'UniformOutput', false );
      %If any row contains an entry that cannot be converted to a number then return the whole
      %column as a char array
      if ~any( cellfun( @isnan, colData ) )
        if ~p.noHeader
          data(:,ii)= cat( 1, data(1,ii), colData{:} );
        else
          data(:,ii)= colData;
        end
      end
    end
  end
  varargout{1}    = data;
else
  %In column mode get rid of the headers (if present)
  data            = data(startRowIdx:end,colIdx);
  if ~p.outputAsChar
    %If we're not outputting the data as characters then try to convert each column to a number
    for ii = 1 : nHeaders
      colData     = cellfun( @str2double, data(:,ii), 'UniformOutput', false );
      %If any row contains an entry that cannot be converted to a number then return the whole
      %column as a char array
      if ~any( cellfun( @isnan, colData ) )
        data(:,ii)= colData;
      else
        %If any column cannot be converted to a number then we cannot convert the output to an array
        %or matrix i.e. uniform output is not possible
        uniformOutputPossible = false;
      end
    end
  end
  if nargout == nHeaders
    %Loop through each column and convert to matrix if possible
    for ii = 1 : nHeaders
      if p.uniformOutput && ~any( cellfun( @ischar, data(:,ii) ) )
        varargout{ii} = cell2mat( data(:,ii) );
      else
        varargout{ii} = data(:,ii);
      end
    end
  else
    %Convert entire table to matrix if possible
    if p.uniformOutput && uniformOutputPossible
      data        =  cell2mat( data );
    end
    varargout{1}  = data;
  end
end
end            
            


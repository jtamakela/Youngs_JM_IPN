%Calculates equilibrium and dynamic Young's moduli for 4 step
%stress-relaxation

% /media/janne/Transcend/Janne/Tyot/Academic/PostDoc/Measurements/IPN/40

clear all
clc; 
close all

%Directory
cd /media/janne/Transcend/Janne/Tyot/Academic/PostDoc/Measurements/IPN/60

% % %% %% %% %% %% %%
Sampling_freq = 10;
% % %% %% %% %% %% %%
r = 0.007/2; %RADIUS
% % %% %% %% %% %% %%


fn = dir('*.CSV'); %Finds .csv files which are then analyzed

no_of_headerlines = 5; %Bose adds 5 lines to csv-file

%h_wait = waitbar(0,'Please wait...');

scrsz = get(0,'ScreenSize');

% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
%i=1;
for i = 1:length(fn); %This defines that whole folder is analyzed.
% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
    
%waitbar(i/length(fn));
% fid = fopen(fn(i).name); %Opens file1 for read access
% %Read the data from the data file, convert, and write to cell array called 'data':        
% data = textscan(fid,'%n%n%n%n%n%n%n','delimiter',',','headerLines',no_of_headerlines);
% fclose(fid);

% OR % % % % % % 
name = fn(i).name;
raw_data = importdata(name);
data = raw_data.data;

F = -data(:,5);
h = -data(:,4);
t = data(:,2);
index = data(:,1);



% figure;
% % plot(t/60, -F);
% % hold on;
% plotyy(t/60, F,t/60, h); %Plotting time in minutes
% 


% Area
A = pi*r^2;


%Stress (Pressure) P=F/A
P = F./A;


figure('Position',[scrsz(3)/2 scrsz(4)/2 scrsz(3)/3 scrsz(4)/3])

subplot(2,1,1)
ax = plotyy(t, P,t, h,'plot'); %Plotting time in seconds

hold on;
title(['Stress-Relaxation measurement, Sample ', fn(i).name],'interpreter','none');

set(get(ax(1),'Ylabel'),'String','Stress (Pa)') 
set(get(ax(2),'Ylabel'),'String','Displacement (mm)') 

xlabel('Time (s)');

%
% P_dyn Stress peaks
% P_eq Stress at equilibrium

max_ind = length(P); %Size of measurement
max_windows = [1 round(max_ind*1/4) round(max_ind*1/2) round(max_ind*3/4) max_ind]; %Windows for peak values

% Peak values and indexes
[P_dyn(1) ind_dyn(1)] = max(P(1:max_windows(2)));
[P_dyn(2) ind_dyn(2)] = max(P(max_windows(2):max_windows(3)));
[P_dyn(3) ind_dyn(3)] = max(P(max_windows(3):max_windows(4)));
[P_dyn(4) ind_dyn(4)] = max(P);
ind_dyn = [ind_dyn(1) ind_dyn(2)+max_windows(2) ind_dyn(3)+max_windows(3) ind_dyn(4)];


% SIMO, NEWPORTILLE VARMAAN SOPII TÄÄ PAREMMIN. ETTII MINIMIT ENNEN
% PIIKKEJÄ. KORVAA ALLAOLEVAT 3-RIVIÄ NÄILLÄ TARVITTAESSA
% Possibility to take equilibrium points from local minimums. Static noise
% makes this less efficient method than latter (Difference ~2%)
% [Low(1) ind_eq(1)] = min(P(1:ind_dyn(1)));
% [Low(2) ind_eq(2)] = min(P(ind_dyn(1):ind_dyn(2)));
% [Low(3) ind_eq(3)] = min(P(ind_dyn(2):ind_dyn(3)));
% [Low(4) ind_eq(4)] = min(P(ind_dyn(3):ind_dyn(4)));
% [Low(5) ind_eq(5)] = min(P(ind_dyn(4):end));
% ind_eq = [ind_eq(1) ind_eq(2)+ind_dyn(1) ind_eq(3)+ind_dyn(2) ind_eq(4)+ind_dyn(3) ind_eq(5)+ind_dyn(4)]; %Indekses where minimums are

plot(t(ind_dyn), P_dyn, 'bo','markersize', 10, 'LineWidth', 3);

% Take the values before straining
%15.15 second strain time (0.05%/0.0033%/s)
ind_eq_end = [round(ind_dyn - 17*Sampling_freq) max_ind]; %step time (~16) + 1 second
ind_eq_Start = ind_eq_end-10*Sampling_freq; %10 window for mean P calculation
ind_eq = [mean([ind_eq_Start; ind_eq_end])];

% %% %% %% 

for k = 1:5;
    P_eq(k) = mean(P(ind_eq_Start(k):ind_eq_end(k))); %Calculate mean P values at equilibrium instead of taking one value
end


plot(t(ind_eq), P_eq, 'ro','markersize', 10, 'LineWidth', 3);


% %%
% Taring and zeroing
% Calculating thickness based on displacement

%Tare 1 F_eq out
P = P-P_eq(1);

for k = 1:5;
    h_eq(k) = mean(h(ind_eq_Start(k)-1000:ind_eq_end(k))); %Calculate mean h from 100 seconds to equilibrium
end

h_zeroed = h_eq-h_eq(1); %Subtracting 0-point before actual measurements


% %% Assumption made -> strain 
strain = [0 0.05 0.1 0.15 0.2];

%Thickness based on the assumption that strains are 0-0.2. 
h_plug = h_zeroed/strain;

% %% 
% Least-squares fitting for Pressure. Young's modulus calculated from that then

S = [strain', ones(length(strain),1)];
E_eq = (S\P_eq');

plot(t(ind_eq), strain.*E_eq(1)+E_eq(2), 'g', 'LineWidth', 2);

% And the same without 1 point %%%%%%%%%%%%%%%%%%%%%%%%

S = [strain(2:end)', ones(length(strain(2:end)),1)];
E_eq_w1 = (S\P_eq(2:end)');

plot(t(ind_eq(2:end)), strain(2:end).*E_eq_w1(1)+E_eq_w1(2), 'r--', 'LineWidth', 2);


% Dynamic Moduli %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P_diff = [P_dyn - P_eq(1:end-1)];
E_dyn = P_diff./0.05; %Strain 5%


%Plotting %% %% %%
subplot(2,1,2)
plot(strain*100, P_eq./1e3, 'ro','markersize', 10, 'LineWidth', 3);
title('Equilibrium stiffnesses')
ylabel('Stress [kPa]');
xlabel('Strain [%]');
hold on;
plot(strain*100, (E_eq(1)./1e3)*strain+(E_eq(2)./1e3), 'g', 'LineWidth', 2);
plot(strain*100, (E_eq_w1(1)./1e3)*strain+(E_eq_w1(2)./1e3), 'r', 'LineWidth', 2);
legend('equilibrium stiffness', 'LS-Fit', 'Fit without 1. value', 'Location', 'NW')


% %% 
% Displaying values
disp(' ')
disp(['Sample ', fn(i).name]);
disp(['Calculated thickness = ', num2str(h_plug), ' mm'])
disp(' ')
disp(['Equilibrium Modulus = ', num2str(E_eq(1)/1e6), ' MPa'])
disp(['Equilibrium Modulus without 1. point = ', num2str(E_eq_w1(1)/1e6), ' MPa'])
disp(' ')
disp(['Dynamic moduli = ', num2str(E_dyn/1e6), ' MPa'])
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('');

%Results %% %% %% %%
Final_name{i,:} = fn(i).name;
Final_thickness(i,:) = h_plug;
Final_E_eq(i,:) = E_eq(1);
Final_E_eq_w1(i,:) = E_eq_w1(1);
Final_E_dyn(i,1:4) = E_dyn;

Final_all = [Final_thickness Final_E_eq Final_E_eq_w1 Final_E_dyn];


%close(h_wait); %Close waitbar


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Steps

%Strain rate is 5um/s

strain_time = h_plug*0.5/0.05;

for i = 1:4;
    step_beginning(i) = max(find(t<t(ind_dyn(i))-strain_time));
end

steptimes = [t(step_beginning); t(end)];

figure; 
for i = 1:4;
subplot(2,2,i);
plot(t,h);
hold on;
plot(t(ind_dyn(i)), h(ind_dyn(i)), 'ro','markersize', 10, 'LineWidth', 3);

plot(t(step_beginning(i)), h(step_beginning(i)), 'go','markersize', 10, 'LineWidth', 3);
xlim([steptimes(i)-5 t(ind_dyn(i))+5]);
title(['step ', num2str(i)]);
end


%plot(t(ind_dyn), P_dyn, 'bo','markersize', 10, 'LineWidth', 3);

disp('');
disp(['Displacements start at = ', num2str(steptimes')])
disp(['Step time = ', num2str(strain_time)])
disp(['Step size = ', num2str(0.05*h_plug), 'mm']) 
disp('');

StepsForInp = [steptimes(1) steptimes(2)-steptimes(1)-strain_time...
    steptimes(3)-steptimes(2)-strain_time steptimes(4)-steptimes(3)-strain_time...
    steptimes(5)-steptimes(4)-strain_time];



%run('/media/janne/Transcend/Janne/Tyot/Academic/PostDoc/Measurements/IPN/Models/Mallit/malliuusiks.m');
% This was done for new IPN analysis
run('/media/janne/Transcend/Janne/Tyot/Academic/PostDoc/Measurements/IPN/Models/Mallit/SameTime_malliuusiks.m');


% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
end
% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
%% 

















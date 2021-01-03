function PLOT_T_PROFILE()
% Usage: PLOT_T_PROFILE()
%
% Purpose: 
%   Load and plot Temperature profile.
%
% Input:
%
% Output:
%
% JMT Nov 2018

pdir = pwd;
cd('..');
addpath([pwd '/mfiles_M3TET']);
addpath([pwd '/mfiles_SPH']);
cd(pdir);

% Specify output directory for the plume model(folder where data is located)
outdir_csv = '/Users/jorge/Tests/SPH/TESTS_SOUTH_ATLANTIC_GPLATES/L1F15/_Subdomains/Ridge_profiles/Temp_130km_depth/_CSV_files';
outdir_csv = [outdir_csv filesep];
outdir_mat = '/Users/jorge/Tests/SPH/TESTS_SOUTH_ATLANTIC_GPLATES/L1F15/_Subdomains/Ridge_profiles';
outdir_mat = [outdir_mat filesep];

cd(outdir_csv);

T.Age_00Ma.A_Ap = csvread('T_profile_1_A_Ap_0Ma.csv');
T.Age_00Ma.B_Bp = csvread('T_profile_2_B_Bp_0Ma.csv');

T.Age_02Ma.A_Ap = csvread('T_profile_1_A_Ap_2Ma.csv');
T.Age_02Ma.B_Bp = csvread('T_profile_2_B_Bp_2Ma.csv');

T.Age_04Ma.A_Ap = csvread('T_profile_1_A_Ap_4Ma.csv');
T.Age_04Ma.B_Bp = csvread('T_profile_2_B_Bp_4Ma.csv');

T.Age_06Ma.A_Ap = csvread('T_profile_1_A_Ap_6Ma.csv');
T.Age_06Ma.B_Bp = csvread('T_profile_2_B_Bp_6Ma.csv');

T.Age_08Ma.A_Ap = csvread('T_profile_1_A_Ap_8Ma.csv');
T.Age_08Ma.B_Bp = csvread('T_profile_2_B_Bp_8Ma.csv');

T.Age_10Ma.A_Ap = csvread('T_profile_1_A_Ap_10Ma.csv');
T.Age_10Ma.B_Bp = csvread('T_profile_2_B_Bp_10Ma.csv');

T.Age_12Ma.A_Ap = csvread('T_profile_1_A_Ap_12Ma.csv');
T.Age_12Ma.B_Bp = csvread('T_profile_2_B_Bp_12Ma.csv');

T.Age_14Ma.A_Ap = csvread('T_profile_1_A_Ap_14Ma.csv');
T.Age_14Ma.B_Bp = csvread('T_profile_2_B_Bp_14Ma.csv');

T.Age_16Ma.A_Ap = csvread('T_profile_1_A_Ap_16Ma.csv');
T.Age_16Ma.B_Bp = csvread('T_profile_2_B_Bp_16Ma.csv');

T.Age_18Ma.A_Ap = csvread('T_profile_1_A_Ap_18Ma.csv');
T.Age_18Ma.B_Bp = csvread('T_profile_2_B_Bp_18Ma.csv');

T.Age_20Ma.A_Ap = csvread('T_profile_1_A_Ap_20Ma.csv');
T.Age_20Ma.B_Bp = csvread('T_profile_2_B_Bp_20Ma.csv');

% Save data
save([outdir_mat 'T_profiles'],'T');

% Plots
field_age     = fieldnames(T);
field_profile = fieldnames(T.Age_00Ma);

j = 0;
k = 1;
for i = 1:11
    
    figure(k); clf
    hold on
    
    subplot(2,1,1)
    plot(T.(field_age{i}).(field_profile{1})(:,2),T.(field_age{i}).(field_profile{1})(:,1),'r-','LineWidth',1.5)
    axis([0 1100 750 1450])
    xlabel('distance (km)')
    ylabel('Temperature (°C)')
    title(['A-A'' ', num2str(j),' Myr'])
    set(gca,'FontSize',16)
    
    subplot(2,1,2)
    plot(T.(field_age{i}).(field_profile{2})(:,2),T.(field_age{i}).(field_profile{2})(:,1),'r-','LineWidth',1.5)
    axis([0 1100 750 1450])
    xlabel('distance (km)')
    ylabel('Temperature (°C)')
    title(['B-B'' ', num2str(j),' Myr'])
    set(gca,'FontSize',16)
    
    j = j+2;
    k = k+1;
end

figure(12);clf
time = [0 2 4 6 8 10 12 14 16 18 20];
for i = 1:11
    subplot(2,1,1)
    hold on
    x  = T.(field_age{i}).(field_profile{1})(:,2)';
    y  = time(i)*ones(size(T.(field_age{i}).(field_profile{1})(:,2),1),1)';
    z  = T.(field_age{i}).(field_profile{1})(:,1)';
    xx =[x;x];
    yy =[y;y];
    zz =[z;z];
    z0 =zeros(size(xx));
    surf(xx,yy,z0,zz,'EdgeColor','interp','linew',15) %// color binded to "y" values
    colormap('jet')
    colorbar
    caxis([750 1450])
    view(2) %// view(0,90)
    axis([0 1100 -0.5 20.5])
    yticks([0 2 4 6 8 10 12 14 16 18 20])
    xlabel('distance (km)')
    ylabel('time (Myr)')
    title('A-A''')
    set(gca,'FontSize',16)
    
    subplot(2,1,2)
    hold on
    x  = T.(field_age{i}).(field_profile{2})(:,2)';
    y  = time(i)*ones(size(T.(field_age{i}).(field_profile{2})(:,2),1),1)';
    z  = T.(field_age{i}).(field_profile{2})(:,1)';
    xx =[x;x];
    yy =[y;y];
    zz =[z;z];
    z0 =zeros(size(xx));
    surf(xx,yy,z0,zz,'EdgeColor','interp','linew',15) %// color binded to "y" values
    colormap('jet')
    colorbar
    caxis([750 1450])
    view(2) %// view(0,90)
    axis([0 1100 -0.5 20.5])
    yticks([0 2 4 6 8 10 12 14 16 18 20])
    xlabel('distance (km)')
    ylabel('time (Myr)')
    title('B-B''')
    set(gca,'FontSize',16)
end

end % END OF FUNCTION PLOT_T_PROFILE
function READ_LITH_DEPTH_CSV_DATA_AND_SAVE_AS_MAT_DATA()
% Usage: READ_LITH_DEPTH_CSV_DATA_AND_SAVE_AS_MAT_DATA()
%
% Purpose: 
%   Read the lithospere depth data that is stored as .csv and save it as
%   .mat data.
%
% Input:
%
% Output:
%
% JMT Oct 2018

pdir = pwd;
cd('..');
addpath([pwd '/mfiles_M3TET']);
addpath([pwd '/mfiles_SPH']);
cd(pdir);

% Specify output directory for the plume model(folder where data is located)
outdir_csv = '/Users/jorge/Tests/SPH/TESTS_SOUTH_ATLANTIC_GPLATES/NoPlume/_Subdomains/Ridge_profiles/Stretching_factor_Beta/_CSV_files';
outdir_csv = [outdir_csv filesep];
outdir_mat = '/Users/jorge/Tests/SPH/TESTS_SOUTH_ATLANTIC_GPLATES/NoPlume/_Subdomains/Ridge_profiles';
outdir_mat = [outdir_mat filesep];

cd(outdir_csv);

% 2 Ma
Depth.Age_1_9Ma.A_Ap = csvread('Depth_1_9Ma_A_Ap.csv');
Depth.Age_1_9Ma.A_B  = csvread('Depth_1_9Ma_A_B.csv');
Depth.Age_1_9Ma.B_Bp = csvread('Depth_1_9Ma_B_Bp.csv');

Depth.Age_2_0Ma.A_Ap = csvread('Depth_2Ma_A_Ap.csv');
Depth.Age_2_0Ma.A_B  = csvread('Depth_2Ma_A_B.csv');
Depth.Age_2_0Ma.B_Bp = csvread('Depth_2Ma_B_Bp.csv');

Depth.Age_2_1Ma.A_Ap = csvread('Depth_2_1Ma_A_Ap.csv');
Depth.Age_2_1Ma.A_B  = csvread('Depth_2_1Ma_A_B.csv');
Depth.Age_2_1Ma.B_Bp = csvread('Depth_2_1Ma_B_Bp.csv');

% 4 Ma
Depth.Age_3_9Ma.A_Ap = csvread('Depth_3_9Ma_A_Ap.csv');
Depth.Age_3_9Ma.A_B  = csvread('Depth_3_9Ma_A_B.csv');
Depth.Age_3_9Ma.B_Bp = csvread('Depth_3_9Ma_B_Bp.csv');

Depth.Age_4_0Ma.A_Ap = csvread('Depth_4Ma_A_Ap.csv');
Depth.Age_4_0Ma.A_B  = csvread('Depth_4Ma_A_B.csv');
Depth.Age_4_0Ma.B_Bp = csvread('Depth_4Ma_B_Bp.csv');

Depth.Age_4_1Ma.A_Ap = csvread('Depth_4_1Ma_A_Ap.csv');
Depth.Age_4_1Ma.A_B  = csvread('Depth_4_1Ma_A_B.csv');
Depth.Age_4_1Ma.B_Bp = csvread('Depth_4_1Ma_B_Bp.csv');

% 6 Ma
Depth.Age_5_9Ma.A_Ap = csvread('Depth_5_9Ma_A_Ap.csv');
Depth.Age_5_9Ma.A_B  = csvread('Depth_5_9Ma_A_B.csv');
Depth.Age_5_9Ma.B_Bp = csvread('Depth_5_9Ma_B_Bp.csv');

Depth.Age_6_0Ma.A_Ap = csvread('Depth_6Ma_A_Ap.csv');
Depth.Age_6_0Ma.A_B  = csvread('Depth_6Ma_A_B.csv');
Depth.Age_6_0Ma.B_Bp = csvread('Depth_6Ma_B_Bp.csv');

Depth.Age_6_1Ma.A_Ap = csvread('Depth_6_1Ma_A_Ap.csv');
Depth.Age_6_1Ma.A_B  = csvread('Depth_6_1Ma_A_B.csv');
Depth.Age_6_1Ma.B_Bp = csvread('Depth_6_1Ma_B_Bp.csv');

% 8 Ma
Depth.Age_7_9Ma.A_Ap = csvread('Depth_7_9Ma_A_Ap.csv');
Depth.Age_7_9Ma.A_B  = csvread('Depth_7_9Ma_A_B.csv');
Depth.Age_7_9Ma.B_Bp = csvread('Depth_7_9Ma_B_Bp.csv');

Depth.Age_8_0Ma.A_Ap = csvread('Depth_8Ma_A_Ap.csv');
Depth.Age_8_0Ma.A_B  = csvread('Depth_8Ma_A_B.csv');
Depth.Age_8_0Ma.B_Bp = csvread('Depth_8Ma_B_Bp.csv');

Depth.Age_8_1Ma.A_Ap = csvread('Depth_8_1Ma_A_Ap.csv');
Depth.Age_8_1Ma.A_B  = csvread('Depth_8_1Ma_A_B.csv');
Depth.Age_8_1Ma.B_Bp = csvread('Depth_8_1Ma_B_Bp.csv');

% 10 Ma
Depth.Age_9_9Ma.A_Ap = csvread('Depth_9_9Ma_A_Ap.csv');
Depth.Age_9_9Ma.A_B  = csvread('Depth_9_9Ma_A_B.csv');
Depth.Age_9_9Ma.B_Bp = csvread('Depth_9_9Ma_B_Bp.csv');

Depth.Age_10_0Ma.A_Ap = csvread('Depth_10Ma_A_Ap.csv');
Depth.Age_10_0Ma.A_B  = csvread('Depth_10Ma_A_B.csv');
Depth.Age_10_0Ma.B_Bp = csvread('Depth_10Ma_B_Bp.csv');

Depth.Age_10_1Ma.A_Ap = csvread('Depth_10_1Ma_A_Ap.csv');
Depth.Age_10_1Ma.A_B  = csvread('Depth_10_1Ma_A_B.csv');
Depth.Age_10_1Ma.B_Bp = csvread('Depth_10_1Ma_B_Bp.csv');

% 12 Ma
Depth.Age_11_9Ma.A_Ap = csvread('Depth_11_9Ma_A_Ap.csv');
Depth.Age_11_9Ma.A_B  = csvread('Depth_11_9Ma_A_B.csv');
Depth.Age_11_9Ma.B_Bp = csvread('Depth_11_9Ma_B_Bp.csv');

Depth.Age_12_0Ma.A_Ap = csvread('Depth_12Ma_A_Ap.csv');
Depth.Age_12_0Ma.A_B  = csvread('Depth_12Ma_A_B.csv');
Depth.Age_12_0Ma.B_Bp = csvread('Depth_12Ma_B_Bp.csv');

Depth.Age_12_1Ma.A_Ap = csvread('Depth_12_1Ma_A_Ap.csv');
Depth.Age_12_1Ma.A_B  = csvread('Depth_12_1Ma_A_B.csv');
Depth.Age_12_1Ma.B_Bp = csvread('Depth_12_1Ma_B_Bp.csv');

% 14 Ma
Depth.Age_13_9Ma.A_Ap = csvread('Depth_13_9Ma_A_Ap.csv');
Depth.Age_13_9Ma.A_B  = csvread('Depth_13_9Ma_A_B.csv');
Depth.Age_13_9Ma.B_Bp = csvread('Depth_13_9Ma_B_Bp.csv');

Depth.Age_14_0Ma.A_Ap = csvread('Depth_14Ma_A_Ap.csv');
Depth.Age_14_0Ma.A_B  = csvread('Depth_14Ma_A_B.csv');
Depth.Age_14_0Ma.B_Bp = csvread('Depth_14Ma_B_Bp.csv');

Depth.Age_14_1Ma.A_Ap = csvread('Depth_14_1Ma_A_Ap.csv');
Depth.Age_14_1Ma.A_B  = csvread('Depth_14_1Ma_A_B.csv');
Depth.Age_14_1Ma.B_Bp = csvread('Depth_14_1Ma_B_Bp.csv');

% 16 Ma
Depth.Age_15_9Ma.A_Ap = csvread('Depth_16Ma_A_Ap.csv');
Depth.Age_15_9Ma.A_B  = csvread('Depth_16Ma_A_B.csv');
Depth.Age_15_9Ma.B_Bp = csvread('Depth_16Ma_B_Bp.csv');

Depth.Age_16_0Ma.A_Ap = csvread('Depth_16_1Ma_A_Ap.csv');
Depth.Age_16_0Ma.A_B  = csvread('Depth_16_1Ma_A_B.csv');
Depth.Age_16_0Ma.B_Bp = csvread('Depth_16_1Ma_B_Bp.csv');

Depth.Age_16_1Ma.A_Ap = csvread('Depth_16_1Ma_A_Ap.csv');
Depth.Age_16_1Ma.A_B  = csvread('Depth_16_1Ma_A_B.csv');
Depth.Age_16_1Ma.B_Bp = csvread('Depth_16_1Ma_B_Bp.csv');

% 18 Ma
Depth.Age_17_9Ma.A_Ap = csvread('Depth_17_9Ma_A_Ap.csv');
Depth.Age_17_9Ma.A_B  = csvread('Depth_17_9Ma_A_B.csv');
Depth.Age_17_9Ma.B_Bp = csvread('Depth_17_9Ma_B_Bp.csv');

Depth.Age_18_0Ma.A_Ap = csvread('Depth_18Ma_A_Ap.csv');
Depth.Age_18_0Ma.A_B  = csvread('Depth_18Ma_A_B.csv');
Depth.Age_18_0Ma.B_Bp = csvread('Depth_18Ma_B_Bp.csv');

Depth.Age_18_1Ma.A_Ap = csvread('Depth_18_1Ma_A_Ap.csv');
Depth.Age_18_1Ma.A_B  = csvread('Depth_18_1Ma_A_B.csv');
Depth.Age_18_1Ma.B_Bp = csvread('Depth_18_1Ma_B_Bp.csv');

% 20 Ma
Depth.Age_19_9Ma.A_Ap = csvread('Depth_19_9Ma_A_Ap.csv');
Depth.Age_19_9Ma.A_B  = csvread('Depth_19_9Ma_A_B.csv');
Depth.Age_19_9Ma.B_Bp = csvread('Depth_19_9Ma_B_Bp.csv');

Depth.Age_20_0Ma.A_Ap = csvread('Depth_20Ma_A_Ap.csv');
Depth.Age_20_0Ma.A_B  = csvread('Depth_20Ma_A_B.csv');
Depth.Age_20_0Ma.B_Bp = csvread('Depth_20Ma_B_Bp.csv');

% Make sure in the initial lithospheric thickness at the southern part is ~235 km

% Random numbers (around depth = 235 km) between dist = 1500 km and dist = 2000 km
min_depth_1_9Ma_B_Bp = min(Depth.Age_1_9Ma.B_Bp(:,1));
max_depth_1_9Ma_B_Bp = max(Depth.Age_1_9Ma.B_Bp(:,1));
r_1_9Ma = min_depth_1_9Ma_B_Bp + ((max_depth_1_9Ma_B_Bp-1)-(min_depth_1_9Ma_B_Bp+1)).*rand(53,1);

depth_1_9Ma_A_Ap              = Depth.Age_1_9Ma.A_Ap(:,1);
depth_1_9Ma_A_Ap_new          = depth_1_9Ma_A_Ap;
depth_1_9Ma_A_Ap_new(148:end) = r_1_9Ma;

% Random numbers (linspace) between dist = 1300 km and dist = 1500 km
depth_1_9Ma_A_Ap_new(128:147) = linspace(126,234,20)' -6 + (6+6).*rand(20,1);

% Create the profiles for 2Myr and 2.1Myr
depth_2_0Ma_A_Ap_new = Depth.Age_2_0Ma.A_Ap(:,1);
depth_2_0Ma_A_Ap_new(128:end) = depth_1_9Ma_A_Ap_new(128:end) - rand(73,1);

depth_2_1Ma_A_Ap_new = Depth.Age_2_1Ma.A_Ap(:,1);
depth_2_1Ma_A_Ap_new(128:end) = depth_2_0Ma_A_Ap_new(128:end) - rand(73,1);

% Plot 
dist_1_9Ma_A_Ap = Depth.Age_1_9Ma.A_Ap(:,2);
dist_1_9Ma_A_B  = Depth.Age_1_9Ma.A_B(:,2);
dist_1_9Ma_B_Bp = Depth.Age_1_9Ma.B_Bp(:,2);

depth_1_9Ma_A_B  = Depth.Age_1_9Ma.A_B(:,1);
depth_1_9Ma_B_Bp = Depth.Age_1_9Ma.B_Bp(:,1);

figure(10);clf
plot(dist_1_9Ma_A_Ap,depth_1_9Ma_A_Ap,'b','LineWidth',1)
hold on
plot(dist_1_9Ma_A_Ap,depth_1_9Ma_A_Ap_new,'r-','LineWidth',1)
% plot(dist_1_9Ma_A_Ap,depth_2_0Ma_A_Ap_new,'g:','LineWidth',2)
% plot(dist_1_9Ma_A_Ap,depth_2_1Ma_A_Ap_new,'m-.','LineWidth',2)
plot(dist_1_9Ma_A_B,depth_1_9Ma_A_B,'b','LineWidth',2)
plot(dist_1_9Ma_B_Bp,depth_1_9Ma_B_Bp,'b','LineWidth',2)
plot(dist_1_9Ma_A_B,depth_1_9Ma_A_B,'r--','LineWidth',2)
plot(dist_1_9Ma_B_Bp,depth_1_9Ma_B_Bp,'r--','LineWidth',2)
grid on
axis([-1100 2100 0 250])
xlabel('distance (km)')
ylabel('lithosphere thickness (km)')
title('Initial lithosphere thickness (2 Myr)')
set(gca,'FontSize',16)
legend('Profile from TecPlot','After editing manually the data','Location','SouthEast')

% Overwrite new profiles
Depth.Age_1_9Ma.A_Ap(:,1) = depth_1_9Ma_A_Ap_new;
Depth.Age_2_0Ma.A_Ap(:,1) = depth_2_0Ma_A_Ap_new;
Depth.Age_2_1Ma.A_Ap(:,1) = depth_2_1Ma_A_Ap_new;

% Save data
save([outdir_mat 'Rigde_profiles'],'Depth');

end % END OF FUNCTION READ_LITH_DEPTH_CSV_DATA_AND_SAVE_AS_MAT_DATA
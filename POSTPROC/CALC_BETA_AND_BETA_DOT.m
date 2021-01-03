function CALC_BETA_AND_BETA_DOT()
% Usage: CALC_BETA_AND_BETA_DOT()
%
% Purpose: 
%   Compute beta and beta_dot from the lithospere thickness data.
%
%   	                   inital lithospere thickness (2 Myr)
%   beta (at 6 Myr)    = --------------------------------------- - 1
%                         current lithosphere thickness (X Myr)
%
%
%                         (beta(6Myr) - beta(5.9Myr)   (beta(6.1Myr) - beta(6Myr)
%                         -------------------------- + --------------------------
%                                    0.1                            0.1
%   beta_dot(at 6 Myr) = ---------------------------------------------------------
%                                                    2
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

outdir_mat = '/Users/jorge/Tests/SPH/TESTS_SOUTH_ATLANTIC_GPLATES/NoPlume/_Subdomains/Ridge_profiles';
outdir_mat = [outdir_mat filesep];

cd(outdir_mat);
load('Rigde_profiles_No_Plume')

init_depth.A_Ap = Depth.Age_2_0Ma.A_Ap(:,1);
init_depth.A_B  = Depth.Age_2_0Ma.A_B(:,1);
init_depth.B_Bp = Depth.Age_2_0Ma.B_Bp(:,1);

field_age     = fieldnames(Depth);
field_profile = fieldnames(Depth.Age_1_9Ma);
for i = 1:numel(field_age)
    for j = 1:numel(field_profile)
        beta.(field_age{i}).(field_profile{j}) = (init_depth.(field_profile{j})./Depth.(field_age{i}).(field_profile{j})(:,1)) - 1;
        
    end
end

for i = 2:3:29
    for j = 1:numel(field_profile)
        if i == 29
            beta_dot.(field_age{i}).(field_profile{j}) = ...
            ( beta.(field_age{i}).(field_profile{j}) - beta.(field_age{i-1}).(field_profile{j})) / 0.1;
        else
            beta_dot.(field_age{i}).(field_profile{j}) = ...
            ( ( beta.(field_age{i}).(field_profile{j}) - beta.(field_age{i-1}).(field_profile{j}) ) / 0.1 + ...
              ( beta.(field_age{i+1}).(field_profile{j}) - beta.(field_age{i}).(field_profile{j}) ) / 0.1 ) / 2;
        end
    end
end

% Plots
j = 0;
for i = 2:3:29
    j = j+2;
    figure(j); clf
    hold on
    
    yyaxis left
    plot(Depth.(field_age{i}).(field_profile{1})(:,2),Depth.(field_age{i}).(field_profile{1})(:,1),'-','LineWidth',1.5)
    plot(Depth.(field_age{i}).(field_profile{2})(:,2),Depth.(field_age{i}).(field_profile{2})(:,1),'-','LineWidth',1.5)
    plot(Depth.(field_age{i}).(field_profile{3})(:,2),Depth.(field_age{i}).(field_profile{3})(:,1),'-','LineWidth',1.5)
    axis([-1100 2100 0 250])
    xlabel('distance (km)')
    ylabel('lithosphere thickness (km)')
    title([num2str(j),' Myr'])
    set(gca,'FontSize',16)
    
    if i == 17 || i == 20 % 12 Myr and 14 Myr
        yyaxis right
        plot(Depth.(field_age{i}).(field_profile{2})(:,2),beta_dot.(field_age{i}).(field_profile{2}),'-','LineWidth',1.5)
        plot(Depth.(field_age{i}).(field_profile{3})(:,2),beta_dot.(field_age{i}).(field_profile{3}),'-','LineWidth',1.5)
        axis([-1100 2100 0 3])
        ylabel('stretching rate (beta dot) (1/Myr)')
    elseif i == 23 % 16 Myr
        yyaxis right
        plot(Depth.(field_age{i}).(field_profile{2})(:,2),beta_dot.(field_age{i}).(field_profile{2}),'-','LineWidth',1.5)
        axis([-1100 2100 0 3])
        ylabel('stretching rate (beta dot) (1/Myr)')
    elseif i == 26 || i == 29 % 18 Myr and 20 Myr
        yyaxis right
        axis([-1100 2100 0 3])
        ylabel('stretching rate (beta dot) (1/Myr)')
    else % 2 Myr, 4 Myr, 6 Myr, 8 Myr and 10 Myr
        yyaxis right
        plot(Depth.(field_age{i}).(field_profile{1})(:,2),beta_dot.(field_age{i}).(field_profile{1}),'-','LineWidth',1.5)
        plot(Depth.(field_age{i}).(field_profile{2})(:,2),beta_dot.(field_age{i}).(field_profile{2}),'-','LineWidth',1.5)
        plot(Depth.(field_age{i}).(field_profile{3})(:,2),beta_dot.(field_age{i}).(field_profile{3}),'-','LineWidth',1.5)
        axis([-1100 2100 0 3])
        ylabel('stretching rate (beta dot) (1/Myr)')
    end
    
    ax = gca;
    ax.XMinorGrid = 'on';
    ax.GridLineStyle = ':';
    ax.XGrid = 'on';
    ax.GridLineStyle = '--';
end

j = 20;
for i = 2:3:29
    j = j+2;
    figure(j); clf
    hold on

    yyaxis left
    plot(Depth.(field_age{i}).(field_profile{1})(:,2),Depth.(field_age{i}).(field_profile{1})(:,1),'-','LineWidth',1.5)
    plot(Depth.(field_age{i}).(field_profile{2})(:,2),Depth.(field_age{i}).(field_profile{2})(:,1),'-','LineWidth',1.5)
    plot(Depth.(field_age{i}).(field_profile{3})(:,2),Depth.(field_age{i}).(field_profile{3})(:,1),'-','LineWidth',1.5)
    axis([-1100 2100 0 250])
    % grid minor
    xlabel('distance (km)')
    ylabel('lithosphere thickness (km)')
    title([num2str(j-20),' Myr'])
    set(gca,'FontSize',16)
    
    if i == 17 || i == 20 || i == 23 % 12 Myr, 14 Myr and 16 Myr
        yyaxis right
        plot(Depth.(field_age{i}).(field_profile{2})(:,2),beta.(field_age{i}).(field_profile{2}),'-','LineWidth',1.5)
        plot(Depth.(field_age{i}).(field_profile{3})(:,2),beta.(field_age{i}).(field_profile{3}),'-','LineWidth',1.5)
        axis([-1100 2100 0 10])
        ylabel('stretching factor (beta)')
    elseif i == 26 || i == 29 % 18 Myr and 20 Myr
        yyaxis right
        axis([-1100 2100 0 10])
        ylabel('stretching factor (beta)')
    else % 2 Myr, 4 Myr, 6 Myr, 8 Myr and 10 Myr
        yyaxis right
        plot(Depth.(field_age{i}).(field_profile{1})(:,2),beta.(field_age{i}).(field_profile{1}),'-','LineWidth',1.5)
        plot(Depth.(field_age{i}).(field_profile{2})(:,2),beta.(field_age{i}).(field_profile{2}),'-','LineWidth',1.5)
        plot(Depth.(field_age{i}).(field_profile{3})(:,2),beta.(field_age{i}).(field_profile{3}),'-','LineWidth',1.5)
        axis([-1100 2100 0 10])
        ylabel('stretching factor (beta)')
    end
    
    ax = gca;
    ax.XMinorGrid = 'on';
    ax.GridLineStyle = ':';
    ax.XGrid = 'on';
    ax.GridLineStyle = '--';
end

end % END FUNCTION CALC_BETA_AND_BETA_DOT
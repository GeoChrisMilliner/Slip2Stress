
clc
clear
close all 

addpath(genpath('./Subroutines2'))

% Author: C. Milliner, Caltech (2022)

% Use of this code should use the following reference:
% Milliner et al. (2022)
% Fault Friction Derived from Fault Bend Influence on Coseismic Slip During
% the 2019 Ridgecrest Earthquake, JGR

% For copyright details see README 

%% Parameters

Num_Zones = 3;% number of stress domains

epsilon_damping_strengthfactor = 0.5;% (default use 0.5) regul. for damping

Write_out_stress_result = 1;

Damping_Inv_method_L1_1_L2_2_conjgrddesc_3 = 3; % L2 inversion produces better R ratio when inverting Hardebeck slip vecs 

Turn_Damping_off = 0;

Num_Bootstrap_iterations = 1000;% number of bootstrap samples to run (usually run 4000)

Output_dir = 'Ridgecrest_case';

%% Read in slip data:

Data_filename_path_4_inv = load('./Input_slip_data/Slip_data_Milliner.txt');

rake_all_zones = Data_filename_path_4_inv(:,3);
strike_all_zones = Data_filename_path_4_inv(:,4);
dip_all_zones = Data_filename_path_4_inv(:,5);
UTMx_slip_vec = Data_filename_path_4_inv(:,1);
UTMy_slip_vec =Data_filename_path_4_inv(:,2);

Hardebeck_dir = [];
disp('----Inverting my slip vecs ----')


%% Ouput filenames 

Filename_save = strcat('./Stress_Inv_Output_Damp/',Output_dir,'/Stress_Inv_output_NumZones_',num2str(Num_Zones),'.mat');
P_T_plot = strcat('./Stress_Figures_Damp/',Output_dir,'/Stress_inv_',num2str(Num_Zones));
stress_Uncert_plot = strcat('./Stress_Figures_Damp/',Output_dir,'/Stress_inv_Uncert_',num2str(Num_Zones));
R_Plot_histogram_figure_name = strcat('./Stress_Figures_Damp/',Output_dir,'/Stress_R_histo_',num2str(Num_Zones));

% make the dir if it doesnt exist
if isempty(dir(strcat('./Stress_Inv_Output_Damp/',Output_dir,'/')))
    mkdir(strcat('./Stress_Inv_Output_Damp/',Output_dir,'/'))
    mkdir(strcat('./Stress_Figures_Damp/',Output_dir,'/'))
end


%% Create stress domains

x_foreshock = [455155.317967936;453399.955190871;452476.080045048;452014.142472136;451275.042355477;446840.441655524;442221.065926406;442221.065926406;449612.067092995;466334.207232401;466241.819717818;460329.018784548;456356.355657507;454878.155424189];
y_foreshock = [3953676.66406516;3951921.30128809;3950720.26359852;3949242.06336520;3945916.11284024;3942312.99977153;3938894.66173198;3936215.42380909;3934829.61109036;3934921.99860494;3939264.21179031;3947394.31307356;3952106.07631726;3952752.78891933];

Indx_north = find(UTMy_slip_vec >= 3961051);% isolate just northern region
Indx_central = find(UTMy_slip_vec < 3961051 & UTMy_slip_vec > 3.9488e6 & UTMx_slip_vec < 4.5196e5);

Indx_south_logical = inpolygon(UTMx_slip_vec,UTMy_slip_vec,x_foreshock,y_foreshock);
Indx_south = find(Indx_south_logical);

Indx_all{1} = Indx_north;
Indx_all{2} = Indx_central;
Indx_all{3} = Indx_south;

% north zone
strike_all_zones_all{1} = strike_all_zones(Indx_north);
dip_all_zones_all{1} = dip_all_zones(Indx_north);
rake_all_zones_all{1} = rake_all_zones(Indx_north);
% south zone
strike_all_zones_all{2} = strike_all_zones(Indx_central);
dip_all_zones_all{2} = dip_all_zones(Indx_central);
rake_all_zones_all{2} = rake_all_zones(Indx_central);
% foreshock zone
strike_all_zones_all{3} = strike_all_zones(Indx_south);
dip_all_zones_all{3} = dip_all_zones(Indx_south);
rake_all_zones_all{3} = rake_all_zones(Indx_south);

figure
scatter(UTMx_slip_vec(Indx_central),UTMy_slip_vec(Indx_central),40,ones(size(Indx_central))*2,'d','filled')
hold on
scatter(UTMx_slip_vec(Indx_north),UTMy_slip_vec(Indx_north),40,ones(size(Indx_north)),'o','filled')
scatter(UTMx_slip_vec(Indx_south),UTMy_slip_vec(Indx_south),40,ones(size(find(Indx_south)))*3,'s','filled')
% plot(x_foreshock,y_foreshock,'r')

legend('Zone 1','Zone 2','Zone 3')

axis equal
title('Stress Domains')

Zone_name_str{1} = 'North Zone';
Zone_name_str{2} = 'Central Zone';
Zone_name_str{3} = 'South Zone';

%% Damping matrix 

I_matrix = eye(5);
D = [I_matrix,I_matrix*-1,zeros(size(I_matrix));...% connection between north and central zone
    zeros(size(I_matrix)),I_matrix,I_matrix*-1];


if Turn_Damping_off
    epsilon_damping_strengthfactor = 0;
    disp('---- WARNING: Turned Damping off -----')
end

Num_stress_zones = length(Indx_all);

figure
for ii = 1:Num_stress_zones
    subplot(1,3,1)
    scatter(UTMx_slip_vec(Indx_all{ii}),UTMy_slip_vec(Indx_all{ii}),20,strike_all_zones(Indx_all{ii})','o','filled')
    colorbar
    axis equal
    title('Strike')
    hold on 
    
    subplot(1,3,2)
    scatter(UTMx_slip_vec(Indx_all{ii}),UTMy_slip_vec(Indx_all{ii}),20,dip_all_zones(Indx_all{ii})','o','filled')
    colorbar
    axis equal
    title('Dip')
    hold on 
    
    
    subplot(1,3,3)
    scatter(UTMx_slip_vec(Indx_all{ii}),UTMy_slip_vec(Indx_all{ii}),20,rake_all_zones(Indx_all{ii})','o','filled')
    axis equal
    colorbar
    title('rake')
    hold on 
    
end

%% Run Div. criteria

[Slip_vecs_div_deg__RMS_all,Angles_in_zone_ii] = Est_slip_vec_diversity(dip_all_zones_all,strike_all_zones_all,Zone_name_str);


%% Invert
Assume_one_Principal_is_vert = 0;

if Assume_one_Principal_is_vert~=0
    error('The variable Assume_one_Principal_is_vert must = 0 as havent implemented other options yet')
end


[Deviatoric_Stress_tensor_for_zz_zone,Obs_Pred_slip_components,L2_misfit_mean,POVR_mean,Soln_roughness_mean,Model_length_mean,...
    Deviatoric_Stress_tensor_for_zz_zone_pcg] ...
    = linear_stress_inversion_CM_damp_conjgrad(strike_all_zones_all,dip_all_zones_all,rake_all_zones_all,Assume_one_Principal_is_vert,...
    D,epsilon_damping_strengthfactor,Damping_Inv_method_L1_1_L2_2_conjgrddesc_3);
    

%--------------------------------------------------------------------------
%% Perform bootstrap by random replacement
%--------------------------------------------------------------------------

h = waitbar(0,'Running Bootstrap. Please wait...');
    

for ii = 1:Num_Bootstrap_iterations
    
    Perc_done = strcat('Running Bootstrap, calculated: ',num2str(round(ii/Num_Bootstrap_iterations*10)/10*100),'%');
    waitbar(ii/Num_Bootstrap_iterations,h,Perc_done)
    
    % bootstrap each zone
    for zz = 1:Num_stress_zones
        strike = strike_all_zones_all{zz};
        dip = dip_all_zones_all{zz};
        rake = rake_all_zones_all{zz};
        
        %% Perform bootstrap
        
        rand_bootstrap_indx = randi(length(strike),length(strike),1);% draw random  integers from a uniform distribution
        % now create the resampled dataset.
        strike_boot{zz} = strike(rand_bootstrap_indx);
        dip_boot{zz} = dip(rand_bootstrap_indx);
        rake_boot{zz} = rake(rand_bootstrap_indx);
        %----------------------------------------
    end
    
    
    [Stress_tensor_for_zz_zone_boootstrap,Obs_Pred_slip_components_USE_Pred,L2_misfit_mean,POVR_mean,Soln_roughness_mean,Model_length_mean,...
        Stress_tensor_for_zz_zone_boootstrap_pcg] ...
        = linear_stress_inversion_CM_damp_conjgrad(strike_boot,dip_boot,rake_boot,Assume_one_Principal_is_vert,D,epsilon_damping_strengthfactor,Damping_Inv_method_L1_1_L2_2_conjgrddesc_3);
    
    Stress_tensor_bootstrap_95CI{ii,:} = Stress_tensor_for_zz_zone_boootstrap;
    
end
close(h)


%% save result
if Write_out_stress_result
        
        Stress_tensor_for_zz_zone = Deviatoric_Stress_tensor_for_zz_zone;        
                        
        Inversion_parameters.damping = epsilon_damping_strengthfactor;
        Inversion_parameters.Zones = Num_stress_zones;

        Inversion_parameters.Num_Bootstrap_resampling = Num_Bootstrap_iterations;
        Inversion_parameters.Inv_method_L1_1_L2_2 = Damping_Inv_method_L1_1_L2_2_conjgrddesc_3;
        Inversion_parameters.Turn_Damping_off = Turn_Damping_off;
        
        Slip_vec_data.strike_all_zones_all = strike_all_zones_all;
        Slip_vec_data.dip_all_zones_all = dip_all_zones_all;
        Slip_vec_data.rake_all_zones_all = rake_all_zones_all;
        Slip_vec_data.UTMx_here = UTMx_slip_vec;
        Slip_vec_data.UTMy_here = UTMy_slip_vec;
        Slip_vec_data.Zone_name_str = Zone_name_str;
        Slip_vec_data.Indx_all = Indx_all;
        
        save(Filename_save,'Obs_Pred_slip_components','Stress_tensor_for_zz_zone',...
            'Stress_tensor_bootstrap_95CI','Indx_all',...
            'Inversion_parameters','Slip_vec_data');

end
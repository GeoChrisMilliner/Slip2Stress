function [Slip_vecs_div_deg__RMS_all,Angles_in_zone_ii] = Est_slip_vec_diversity(dip_all_zones_all,strike_all_zones_all,Zone_name_str)

Min_OK_RMS_diversity = 25;
    
Num_stress_zones = length(dip_all_zones_all);

figure
for ii =1:Num_stress_zones
    
    dips_in_zone_ii = dip_all_zones_all{ii};
    strikes_in_zone_ii = strike_all_zones_all{ii};

    Indx_adjust_just_4_this = find(strikes_in_zone_ii<180);
    dips_in_zone_ii(Indx_adjust_just_4_this) = 180-dips_in_zone_ii(Indx_adjust_just_4_this);% <-- make it negative.
    strikes_in_zone_ii(Indx_adjust_just_4_this) = strikes_in_zone_ii(Indx_adjust_just_4_this)+180;
    
    Avg_dip = nanmedian(dips_in_zone_ii);
    Avg_strike = nanmedian(strikes_in_zone_ii);
    
    Num_data_points_in_zone_ii = length(strikes_in_zone_ii);
    
    
    subplot(Num_stress_zones,2,ii*2-1)
    histogram(dip_all_zones_all{ii},8)% before correction
    hold on 
    histogram(dips_in_zone_ii,8)
    xline(Avg_dip,'r--')
    title(strcat('Dip Zone',Zone_name_str{ii}));% after correction
    if ii == 1
    legend('Before correction','after','location','best')
    end
    
    
    subplot(Num_stress_zones,2,ii*2)
    histogram(strike_all_zones_all{ii},8)
    hold on 
    histogram(strikes_in_zone_ii,8)
    xline(Avg_strike,'r--')
    title(strcat('Strike Zone: ',Zone_name_str{ii}));
    
    for jj = 1:Num_data_points_in_zone_ii
        jth_strikes_in_zone_ii = strikes_in_zone_ii(jj);
        jth_dip_in_zone_ii = dips_in_zone_ii(jj);
        
        [Angle_between_deg(jj)] = Angle_btwn_2_fault_planes(Avg_strike,Avg_dip,jth_strikes_in_zone_ii,jth_dip_in_zone_ii);
    
    end
    
    Angles_in_zone_ii{ii} = Angle_between_deg;
    % RMS of angle of slip vectors from the mean slip vector of the group
    Slip_vecs_div_deg__RMS = sqrt(nanmean(Angle_between_deg.^2));
    
    if Slip_vecs_div_deg__RMS >=Min_OK_RMS_diversity
        disp('---- Sufficient diversity of slip planes -----')
    elseif Slip_vecs_div_deg__RMS<Min_OK_RMS_diversity
        disp('*!*!**! WARNING diversity of slip planes is NOT sufficient*!*!**!')
    end


    Slip_vecs_div_deg__RMS_all{ii} = Slip_vecs_div_deg__RMS;
    
    clear Angle_between_deg
end


figure

for ii = 1:Num_stress_zones
    histogram(Angles_in_zone_ii{ii})
    hold on 
end
xlabel('Angle (^o) from average slip vector')
set(gca,'fontsize',16)
title({'Slip vector Diversity','Min. RMS should be = 30^o'})
ylabel('Frequency')


legend(strcat(Zone_name_str{1},' RMS = ',num2str(round(Slip_vecs_div_deg__RMS_all{1})),'^o'),strcat(Zone_name_str{2},' RMS = ',num2str(round(Slip_vecs_div_deg__RMS_all{2})),'^o'),strcat(Zone_name_str{3},' RMS = ',num2str(round(Slip_vecs_div_deg__RMS_all{3})),'^o'))



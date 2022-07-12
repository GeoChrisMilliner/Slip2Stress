function [Angle_between_deg] = Angle_btwn_2_vecs(Vector_1_strike_deg,Vector_1_dip_deg,Vector_1_rake_deg,Vector_2_strike_deg,Vector_2_dip_deg,Vector_2_rake_deg)

% find the angle between two 3D vectors, not fault planes! 

%% input should be:
% - fault strike
% - fault dip (not the vector plunge) 
% - vector rake. 

% N.B. all input should be in deg. 


%% Convert strike, dip rake --> E,N,V Cartesian co-ordinates
% see page 24 of: http://www.geol.tsukuba.ac.jp/~yagi-y/text/Source_meca_v0.7.pdf
% or pg 24 of: Eq_for_rake__pg24.pdf
s1_north =  cos(Vector_1_rake_deg*pi/180).*cos(Vector_1_strike_deg*pi/180) + cos(Vector_1_dip_deg*pi/180).*sin(Vector_1_rake_deg*pi/180).*sin(Vector_1_strike_deg*pi/180);%dNorth
s1_east =  cos(Vector_1_rake_deg*pi/180).*sin(Vector_1_strike_deg*pi/180) - cos(Vector_1_dip_deg*pi/180).*sin(Vector_1_rake_deg*pi/180).*cos(Vector_1_strike_deg*pi/180);%dEast
s1_vert = -sin(Vector_1_rake_deg*pi/180).*sin(Vector_1_dip_deg*pi/180);% vertical component

s2_north =  cos(Vector_2_rake_deg*pi/180).*cos(Vector_2_strike_deg*pi/180) + cos(Vector_2_dip_deg*pi/180).*sin(Vector_2_rake_deg*pi/180).*sin(Vector_2_strike_deg*pi/180);%dNorth
s2_east =  cos(Vector_2_rake_deg*pi/180).*sin(Vector_2_strike_deg*pi/180) - cos(Vector_2_dip_deg*pi/180).*sin(Vector_2_rake_deg*pi/180).*cos(Vector_2_strike_deg*pi/180);%dEast
s2_vert = -sin(Vector_2_rake_deg*pi/180).*sin(Vector_2_dip_deg*pi/180);% vertical component

%% angle between two unit normal vectors. 

for ii = 1:length(s2_north)
    unit_vec1 = [s1_north(ii);s1_east(ii);s1_vert(ii)];
    unit_vec2 = [s2_north(ii);s2_east(ii);s2_vert(ii)];
    
    Angle_between_deg(ii) = acosd(unit_vec1'*unit_vec2);
end


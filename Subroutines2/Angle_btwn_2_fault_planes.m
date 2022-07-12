function [Angle_between_deg,Angle_between_deg_matlab] = Angle_btwn_2_fault_planes(Fault_1_strike_deg,Fault_1_dip_deg,Fault_2_strike_deg,Fault_2_dip_deg)

% - find the angle between two fault planes by calc. the arc cosine of the dot prod.
% between two unit vecotrs normal to each plane. 
% - so essentially finding the angle between the two unit vectors that are normal 
% to the fault surface 

% All input should be in degrees. 

%% Get the normal vectors for each fault plane 
n1a = -sin(Fault_1_dip_deg*pi/180).*sin(Fault_1_strike_deg*pi/180);
n2a =  sin(Fault_1_dip_deg*pi/180).*cos(Fault_1_strike_deg*pi/180);
n3a = -cos(Fault_1_dip_deg*pi/180);


n1b = -sin(Fault_2_dip_deg*pi/180).*sin(Fault_2_strike_deg*pi/180);
n2b =  sin(Fault_2_dip_deg*pi/180).*cos(Fault_2_strike_deg*pi/180);
n3b = -cos(Fault_2_dip_deg*pi/180);

%% angle between two unit normal vectors. 
unit_vec1 = [n1a;n2a;n3a];
unit_vec2 = [n1b;n2b;n3b];

Angle_between_deg = acosd(unit_vec1'*unit_vec2);

% MATLAB
% CosTheta = max(min(dot(unit_vec1,unit_vec2)/(norm(unit_vec1)*norm(unit_vec2)),1),-1);
% Angle_between_deg_matlab = real(acosd(CosTheta));


%*************************************************************************%
%                                                                         %
%  function AZIMUTH_PLUNGE                                                %
%                                                                         %
%  calculation of azimuth and plunge of principal stress axes             %
%                                                                         %
%  input:  stress tensor                                                  %
%  output: azimuth and plunge of principal stress axes                    %
%                                                                         %
%*************************************************************************%
function [direction_sigma_1 direction_sigma_2 direction_sigma_3] = azimuth_plunge(tau)

%--------------------------------------------------------------------------
% eigenvalues and eienvectors of the stress tensor
%--------------------------------------------------------------------------
[vector diag_tensor] = eig(tau);

value = eig(diag_tensor);
[value_sorted,j] = sort(value);

sigma_vector_1 = vector(:,j(1));
sigma_vector_2 = vector(:,j(2));
sigma_vector_3 = vector(:,j(3));

if (sigma_vector_1(3)<0) sigma_vector_1 = -sigma_vector_1; end
if (sigma_vector_2(3)<0) sigma_vector_2 = -sigma_vector_2; end
if (sigma_vector_3(3)<0) sigma_vector_3 = -sigma_vector_3; end

sigma = sort(eig(tau));
shape_ratio = (sigma(1)-sigma(2))/(sigma(1)-sigma(3));

%--------------------------------------------------------------------------
% 1. eigenvector
%--------------------------------------------------------------------------
fi_1 = atan(abs(sigma_vector_1(2)./sigma_vector_1(1)))*180/pi;

if (sigma_vector_1(2)>=0 && sigma_vector_1(1)>=0) azimuth_sigma_1 =     fi_1; end  
if (sigma_vector_1(2)>=0 && sigma_vector_1(1)< 0) azimuth_sigma_1 = 180-fi_1; end  
if (sigma_vector_1(2)< 0 && sigma_vector_1(1)< 0) azimuth_sigma_1 = 180+fi_1; end  
if (sigma_vector_1(2)< 0 && sigma_vector_1(1)>=0) azimuth_sigma_1 = 360-fi_1; end  

theta_sigma_1 = acos(abs(sigma_vector_1(3)))*180/pi;
    
%--------------------------------------------------------------------------
% 2. eigenvector
%--------------------------------------------------------------------------
fi_2 = atan(abs(sigma_vector_2(2)./sigma_vector_2(1)))*180/pi;

if (sigma_vector_2(2)>=0 && sigma_vector_2(1)>=0) azimuth_sigma_2 =     fi_2; end  
if (sigma_vector_2(2)>=0 && sigma_vector_2(1)< 0) azimuth_sigma_2 = 180-fi_2; end  
if (sigma_vector_2(2)< 0 && sigma_vector_2(1)< 0) azimuth_sigma_2 = 180+fi_2; end  
if (sigma_vector_2(2)< 0 && sigma_vector_2(1)>=0) azimuth_sigma_2 = 360-fi_2; end  

theta_sigma_2 = acos(abs(sigma_vector_2(3)))*180/pi;

%--------------------------------------------------------------------------
% 3. eigenvector
%--------------------------------------------------------------------------
fi_3 = atan(abs(sigma_vector_3(2)./sigma_vector_3(1)))*180/pi;

if (sigma_vector_3(2)>=0 && sigma_vector_3(1)>=0) azimuth_sigma_3 =     fi_3; end  
if (sigma_vector_3(2)>=0 && sigma_vector_3(1)< 0) azimuth_sigma_3 = 180-fi_3; end  
if (sigma_vector_3(2)< 0 && sigma_vector_3(1)< 0) azimuth_sigma_3 = 180+fi_3; end  
if (sigma_vector_3(2)< 0 && sigma_vector_3(1)>=0) azimuth_sigma_3 = 360-fi_3; end  

theta_sigma_3 = acos(abs(sigma_vector_3(3)))*180/pi;

%--------------------------------------------------------------------------
% azimuth and plunge of the stress axes
%--------------------------------------------------------------------------
direction_sigma_1 = [azimuth_sigma_1 90-theta_sigma_1];
direction_sigma_2 = [azimuth_sigma_2 90-theta_sigma_2];
direction_sigma_3 = [azimuth_sigma_3 90-theta_sigma_3];

end


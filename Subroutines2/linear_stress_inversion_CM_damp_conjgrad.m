function [stress,Obs_Pred_slip_components,L2_misfit,POVR,Soln_roughness,Model_length,stress_vector_all_other] ...
    = linear_stress_inversion_CM_damp_conjgrad(strike_all,dip_all,rake_all,Assume_one_Principal_is_vert,D,epsilon_damping_strengthfactor,Inv_method_L1_1_L2_2_conj_grad_descent)


%--------------------------------------------------------------------------
%  fault normals and slip directions
%--------------------------------------------------------------------------
% this is in E,N,V cartesian co-ord. system

Num_stress_zones = length(dip_all);

if ~iscell(dip_all) || ~iscell(rake_all) || ~iscell(strike_all)
    error('the strike, dip and rake variables must be cell arrays with each entry the values for each stress zone')
end


for jj = 1:Num_stress_zones
    Num_data_points_each_zone(jj) = length(dip_all{jj});
end

Total_num_vecs = sum(Num_data_points_each_zone);

Num_stress_param = 5;
A_all = zeros(Total_num_vecs*3,Num_stress_zones*Num_stress_param);
d_vector = [];

Prev_row = 0;
Prev_col = 0;
for ii = 1:Num_stress_zones
    
    strike = strike_all{ii};
    dip = dip_all{ii};
    rake = rake_all{ii};
    
    strike = strike(:);
    dip = dip(:);
    rake = rake(:);
    
    % components of the slip vector in: [dNorth (s1),dEast (s2), and vertical (s3)].
    s1 =  cos(rake*pi/180).*cos(strike*pi/180) + cos(dip*pi/180).*sin(rake*pi/180).*sin(strike*pi/180);%dNorth
    s2 =  cos(rake*pi/180).*sin(strike*pi/180) - cos(dip*pi/180).*sin(rake*pi/180).*cos(strike*pi/180);%dEast
    s3 = -sin(rake*pi/180).*sin(dip*pi/180);% vertical component
   
    % normal vector
    n1 = -sin(dip*pi/180).*sin(strike*pi/180);
    n2 =  sin(dip*pi/180).*cos(strike*pi/180);
    n3 = -cos(dip*pi/180);
    
    %--------------------------------------------------------------------------
    % inverted matrix A
    %--------------------------------------------------------------------------
    % matrix coefficients
    
    % Matrix(row,col)
    % row 1
    A_11 = n1-n1.^3+n1.*n3.^2;
    A_12 = n2-2*n2.*n1.^2;
    A_13 = n3-2*n3.*n1.^2;
    A_14 = -n1.*n2.^2+n1.*n3.^2;
    A_15 = -2*n1.*n2.*n3;
    
    % row 2
    A_21 = -n2.*n1.^2+n2.*n3.^2;
    A_22 = n1-2*n1.*n2.^2;
    A_23 = -2*n1.*n2.*n3;
    A_24 = n2-n2.^3+n2.*n3.^2;
    A_25 =  n3-2*n3.*n2.^2;
    
    % row 3
    A_31 =  -n3.*n1.^2 -n3+n3.^3;
    A_32 =  -2*n1.*n2.*n3;
    A_33 =  n1-2*n1.*n3.^2;
    A_34 =  -n2.^2.*n3-n3+n3.^3;
    A_35 =  n2-2*n2.*n3.^2;
    	

    A = [A_11,A_12,A_13,A_14,A_15;...
        A_21,A_22,A_23,A_24,A_25;...
        A_31,A_32,A_33,A_34,A_35];
    
    if Assume_one_Principal_is_vert
        % remove the columns associated with the axis that we assume is
        % vertical.
        A(:,5) = [];
        A(:,3) = [];
        
    end
    
    A_all(Prev_row+1:Prev_row+size(A,1),Prev_col+1:Prev_col+size(A,2)) = A;
    Prev_row = Prev_row + size(A,1);
    Prev_col = Prev_col + size(A,2);
    
    A_store{ii}=A;
    s1_store{ii} = s1;
    s2_store{ii} = s2;
    s3_store{ii} = s3;
    
    d_vector = [d_vector; s1(:); s2(:); s3(:)];
    
end

A_all = sparse(A_all);
D = sparse(D);

%--------------------------------------------------------------------------
%% Inversion
%--------------------------------------------------------------------------
eps = epsilon_damping_strengthfactor;
stress_vector_all_other = [];
switch Inv_method_L1_1_L2_2_conj_grad_descent
    case 1
        Plot_fig = 0;
        Fig_val = 200;
        Criteria_val_stop_L1 = 0.01;
        epsilon = 1e-5;
        [stress_vector_all,~] = L1_norm_CM_Menke_Harde_Damped(A_all,d_vector,Criteria_val_stop_L1,epsilon,Plot_fig,Fig_val,D,eps);
        
    case 2
        
        stress_vector_all = inv(A_all'*A_all  + eps^2*(D'*D))*A_all'*d_vector;% eq. 14 Hardebeck and Michael (2006)
        
    case 3
        
%         A x = b <-- solve for this, so reformulate the normal equations. 

        A_effective = A_all'*A_all  + eps^2*(D'*D);
        b_effective = A_all'*d_vector;
        
        maxit = 60;
        tol = 1e-6;
        
        [stress_vector_all,flag_cgs,relres_cgs,iter_cgs] = cgs(A_effective,b_effective,tol,maxit);

end

row_start = 0;
for ii = 1:Num_stress_zones
    
    stress_vector = stress_vector_all(row_start+1:row_start+5);
    row_start = row_start+5;
    
    s1 = s1_store{ii};
    s2 = s2_store{ii};
    s3 = s3_store{ii};
   
    %% Construct stress tensor
    if Assume_one_Principal_is_vert
        
        sigma_33 = -stress_vector(1)-stress_vector(3);% have a constraint of zero trace of the stress tensor.
        
        stress_tensor = [stress_vector(1),stress_vector(2),0;...
            stress_vector(2),stress_vector(3),0;...
            0,0,sigma_33];
        
        disp(' ----- Assuming one of the Principal stresses are vertical ----- ')
    else
        sigma_33 = -stress_vector(1)-stress_vector(4);% have a constraint of zero trace of the stress tensor.
        
        
        stress_tensor = [stress_vector(1),stress_vector(2),stress_vector(3);...
            stress_vector(2),stress_vector(4),stress_vector(5);...
            stress_vector(3),stress_vector(5),sigma_33];
        
    end
    
    %% Prediction
    
%     sigma  = eig(stress_tensor);
%     stress{ii} = stress_tensor/max(abs(sigma));
    stress{ii} = stress_tensor;
    
%     stress_vector = stress_vector/max(abs(sigma));
    
    % stress_vector_CM = inv(A'*A)*A'*a_vector;
    % pred.
    Prediction_u_slip_components_ith = A_store{ii}*stress_vector;
    Prediction_u_slip_components{ii} = Prediction_u_slip_components_ith;
    
    % Pred. components of slip unit vector in east, north and up co-ord system!
    Pred_s1 = Prediction_u_slip_components_ith(1:length(s1(:)));
    Pred_s2 = Prediction_u_slip_components_ith(length(s1(:))+1:length(s1(:))+length(s1(:)));
    Pred_s3 = Prediction_u_slip_components_ith(length(s1(:))+length(s2(:))+1:end);
    
    Obs_Pred_slip_components{ii}.s1 = s1;% dNorth
    Obs_Pred_slip_components{ii}.s2 = s2;% dEast
    Obs_Pred_slip_components{ii}.s3 = s3;% dV
    Obs_Pred_slip_components{ii}.Pred_s1 = Pred_s1;% dNorth
    Obs_Pred_slip_components{ii}.Pred_s2 = Pred_s2;% dEast
    Obs_Pred_slip_components{ii}.Pred_s3 = Pred_s3;% dV
    
    
    
    
end

%% Misfit:

% L2 norm
Pred = A_all*stress_vector_all;
residual = Pred - d_vector;

% eq 13 of Hardebeck and Michael (2005)
L2_misfit = sqrt(sum(residual.^2));

POVR = (1 - nansum((residual.^2))./nansum(d_vector.^2))*100;

Soln_roughness = stress_vector_all(:)'*stress_vector_all(:);

% eq 12 of Hardebeck and Michael (2005)
% Model_length = norm(D*stress_vector_all);
Model_length = sqrt(sum((D*stress_vector_all).^2));


end



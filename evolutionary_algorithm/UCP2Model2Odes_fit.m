function h = UCP2Model2Odes_fit(para,f_dox)

% The copyright belongs to Kang-Ling Liao, Assistant Professor, 
% Departments of Mathematics and Biological Sciences, University of Manitoba.

h = @fx;
    
lambda_C = para(1);
tildeC = para(2);
eta_T = para(3);
d_C = para(4);
tildelambda_DC = para(5);
K_C = para(6);
lambda_DU = para(7);
K_U = para(8);
d_D = para(9);
tildelambda_M_2C = para(10);
lambda_M_2Q = para(11);
K_Q = para(12);
d_M_2 = para(13);
tildelambda_T = para(14);
K_D = para(15);
barK_M_2 = para(16);
barK_Q = para(17);
lambda_TU = para(18);
d_T = para(19);
lambda_UC = para(20);
lambda_UM_2 = para(21);
K_M_2 = para(22);
d_U = para(23);
rho_P = para(24);
rho_L = para(25);
epsilon = para(26);
lambda_LU = para(27);

    function [yout] = fx(t,y)       
        P = rho_P*y(4);
        L = rho_L*(y(4)+epsilon*y(3)+epsilon*y(1))*(1+lambda_LU*y(5));
        yout = zeros(5,1);      
        yout(1) = lambda_C*y(1)*(1-y(1)/tildeC) - eta_T*y(4)*y(1) - d_C*y(1);
        yout(2) = tildelambda_DC*y(1)/(K_C+y(1)) + lambda_DU*y(5)/(K_U + y(5)) - d_D*y(2);
        yout(3) = tildelambda_M_2C*y(1)/(K_C+y(1)) + lambda_M_2Q*y(3)*(1+P*L/(K_Q + P*L)) - d_M_2*y(3);
        yout(4) = tildelambda_T*y(2)/(K_D+y(2))*1/(1+y(3)/barK_M_2)*1/(1+P*L/barK_Q)+lambda_TU*y(5)/(K_U+y(5)) - d_T*y(4);
        yout(5) = f_dox*(lambda_UC*y(1)/(K_C+y(1)) + lambda_UM_2*y(3)/(K_M_2+y(3)) - d_U*y(5));              
    end
end



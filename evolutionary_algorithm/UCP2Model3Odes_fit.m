function h = UCP2Model3Odes_fit(para)

% The copyright belongs to Kang-Ling Liao, Assistant Professor, 
% Departments of Mathematics and Biological Sciences, University of Manitoba.

h = @fx;

lambda_C = para(1);
tildeC = para(2);
eta_T = para(3);
d_C = para(4);
tildelambda_T = para(5);
K_C = para(6);
barK_Q = para(7);
lambda_TU = para(8);
K_U = para(9);
d_T = para(10);
lambda_UC = para(11);
lambda_UQ = para(12);
K_Q = para(13);
d_U = para(14);
rho_P = para(15);
rho_L = para(16);
epsilon = para(17);
lambda_LU = para(18);

    function [yout] = fx(t,y)       
        P = rho_P*y(2);
        L = rho_L*(y(2)+epsilon*y(1))*(1+lambda_LU*y(3));
        yout = zeros(3,1);      
        yout(1) = lambda_C*y(1)*(1-y(1)/tildeC) - eta_T*y(2)*y(1) - d_C*y(1);
        yout(2) = tildelambda_T*y(1)/(K_C+y(1))*1/(1+P*L/barK_Q) + lambda_TU*yout(3)/(K_U+yout(3)) - d_T*y(2);
        yout(3) = lambda_UC*y(1)/(K_C+y(1)) + lambda_UQ*P*L/(K_Q+P*L) - d_U*y(3);
    end
end



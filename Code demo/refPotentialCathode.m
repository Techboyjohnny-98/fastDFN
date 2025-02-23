%% Reference Potential for Pos. Electrode: Unref(theta_n)
%   Created July 12, 2011 by Scott Moura

function [Uref,varargout] = refPotentialCathode(p,theta)
OCP_Cathode_1 = p.OCP_Cathode(1);
OCP_Cathode_2 = p.OCP_Cathode(2);
OCP_Cathode_3 = p.OCP_Cathode(3);
OCP_Cathode_4 = p.OCP_Cathode(4);
OCP_Cathode_5 = p.OCP_Cathode(5);
OCP_Cathode_6 = p.OCP_Cathode(6);
OCP_Cathode_7 = p.OCP_Cathode(7);
OCP_Cathode_8 = p.OCP_Cathode(8);
OCP_Cathode_9 = p.OCP_Cathode(9);
OCP_Cathode_10 = p.OCP_Cathode(10);
OCP_Cathode_11 = p.OCP_Cathode(11);
OCP_Cathode_12 = p.OCP_Cathode(12);
OCP_Cathode_13 = p.OCP_Cathode(13);
OCP_Cathode_14 = p.OCP_Cathode(14);
OCP_Cathode_15 = p.OCP_Cathode(15);
OCP_Cathode_16 = p.OCP_Cathode(16);
OCP_Cathode_17 = p.OCP_Cathode(17);
OCP_Cathode_18 = p.OCP_Cathode(18);
OCP_Cathode_19 = p.OCP_Cathode(19);

% % Polynomail Fit
% Uref = ppvalFast(p.Uppp,theta);

% DUALFOIL: CoO2 (Cobalt Dioxide) 0.5 < y < 0.99
Uref = OCP_Cathode_1...
 +OCP_Cathode_2*tanh(OCP_Cathode_3-OCP_Cathode_4*theta) ...
 + OCP_Cathode_5*tanh(OCP_Cathode_6-OCP_Cathode_7*theta) ...
 - OCP_Cathode_8*tanh(OCP_Cathode_9-OCP_Cathode_10*theta) ...
 + OCP_Cathode_11*tanh(OCP_Cathode_12-OCP_Cathode_13*theta) ...
 + OCP_Cathode_14*tanh((-theta+OCP_Cathode_15)/OCP_Cathode_16) ...
 - OCP_Cathode_17*tanh((theta-OCP_Cathode_18)/OCP_Cathode_19);

% Gradient of OCP wrt theta
if(nargout == 2)

%     % Polynomail Fit
%     dUref = ppvalFast(p.dUppp,theta);
%     varargout{1} = dUref / p.c_s_p_max;

dUref = OCP_Cathode_2*(-OCP_Cathode_4/p.c_s_p_max)* ...
((1.0./cosh(OCP_Cathode_3-OCP_Cathode_4*theta)).^2) ...
+OCP_Cathode_5*(-OCP_Cathode_7/p.c_s_p_max)*((cosh(OCP_Cathode_6-OCP_Cathode_7*theta)).^(-2)) ...
+OCP_Cathode_8*(19.854/p.c_s_p_max)*((cosh(OCP_Cathode_9-OCP_Cathode_10*theta)).^(-2)) ...
-OCP_Cathode_11*(OCP_Cathode_13/p.c_s_p_max)*((cosh(OCP_Cathode_12-OCP_Cathode_13*theta)).^(-2)) ...
-OCP_Cathode_14/OCP_Cathode_16/p.c_s_p_max*((cosh((-theta+OCP_Cathode_15)/OCP_Cathode_16)).^(-2)) ...
 - OCP_Cathode_17/OCP_Cathode_19/p.c_s_p_max*((cosh((theta-OCP_Cathode_18)/OCP_Cathode_19)).^(-2));
varargout{1} = dUref;

end

% Gradient of OCP wrt temperature
if(nargout >= 3)

    dUdT = 0;
    varargout{2} = dUdT;
end
%% Original
% % % Polynomail Fit
% % Uref = ppvalFast(p.Uppp,theta);
% 
% % DUALFOIL: CoO2 (Cobalt Dioxide) 0.5 < y < 0.99
% Uref = 2.16216+0.07645*tanh(30.834-54.4806*theta) ...
%  + 2.1581*tanh(52.294-50.294*theta) ...
%  - 0.14169*tanh(11.0923-19.8543*theta) ...
%  + 0.2051*tanh(1.4684-5.4888*theta) ...
%  + 0.2531*tanh((-theta+0.56478)/0.1316) ...
%  - 0.02167*tanh((theta-0.525)/0.006);
% 
% % Gradient of OCP wrt theta
% if(nargout == 2)
% 
% %     % Polynomail Fit
% %     dUref = ppvalFast(p.dUppp,theta);
% %     varargout{1} = dUref / p.c_s_p_max;
% 
% dUref = 0.07645*(-54.4806/p.c_s_p_max)* ...
% ((1.0./cosh(30.834-54.4806*theta)).^2) ...
% +2.1581*(-50.294/p.c_s_p_max)*((cosh(52.294-50.294*theta)).^(-2)) ...
% +0.14169*(19.854/p.c_s_p_max)*((cosh(11.0923-19.8543*theta)).^(-2)) ...
% -0.2051*(5.4888/p.c_s_p_max)*((cosh(1.4684-5.4888*theta)).^(-2)) ...
% -0.2531/0.1316/p.c_s_p_max*((cosh((-theta+0.56478)/0.1316)).^(-2)) ...
%  - 0.02167/0.006/p.c_s_p_max*((cosh((theta-0.525)/0.006)).^(-2));
% varargout{1} = dUref;
% 
% end
% 
% % Gradient of OCP wrt temperature
% if(nargout >= 3)
% 
%     dUdT = 0;
%     varargout{2} = dUdT;
% 
% end
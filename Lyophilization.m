clear, close, clc
global k h T_P T_F rho_s W delH_s
T_P = 120+273.15;  % Plate Temp - [K]
T_F = -20+273.15;  % Sublimation Temp - [K]
H = 0.03;   % Bed Height - [m]
W = 0.9;    % Moisture - [kg/kg]
sigma = 4.87e-8;    % Boltzman Coefficient [kcal/m2/h/K4]
rho_s = 1000;   % Solid Density - [kg/m3]
delH_s = 680;   % Sublimation Enthalpy - [kcal/kg]
k = 0.1;        % Solid Conductivity - [kcal/m/h/K]
T_S_max = 50+273.15;   % Maximum Solid Temp. [K]

h = sigma*(T_P^4-T_S_max^4)/(T_P-T_S_max);  % [kcal/m2/h/K]
time_1 = lyophilization_time(H);
disp("The complete lyophilization occurs at x = "+H+" m and requires "+time_1+" hours!")

tspan = 0:0.1:50;
[t,x] = ode23s(@distance_ODE,tspan,0);


for i= 1:length(x)
    T_S(i) = (k*T_F+h*x(i)*(T_P-50))/(k+h*x(i));
end

plot(t,T_S-273.15); hold on; plot(t,T_S_max*ones(1,length(x))-273.15); xlabel('time [hr]'); ylabel('T_{S} [degC]')

function outODE = distance_ODE(t, x)
    global k h T_P T_F rho_s W delH_s
    outODE = (k * h / (k + h * x) * (T_P - T_F) / delH_s) / (rho_s * W);
end

function out = lyophilization_time(x)
global k h T_P T_F rho_s W delH_s
arg_1 = k*x+0.5*h*x^2;
arg_2 = k*h*(T_P - T_F)/(rho_s*W*delH_s);
out = arg_1/arg_2;
end
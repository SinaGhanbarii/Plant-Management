clf, clear, close, clc
global P Gdry A B C rhoL Kpw nAvp h MWwater MWair dHev CpL CpG rhoG muG dp mp0 vg
D = 9;
Tp0=32+273.15;                    % Initial Particle Temperature [K]
Tg0=130+273.15;                    % Initial Gas Temperature [K]
P=1;                        % Pressure [atm]
CpL=1;                      % Specific Heat of Liquid - [kcal/kg/K]
CpG=0.25;                   % Specific Heat of Liquid - [kcal/kg/K]
MWwater=18;                 % Molecular Weight of Water - [kg/kmol]
MWair=29;                   % Molecular Weight of Water - [kg/kmol]
R=0.0821;                   % [m^3*atm/K/kmol]
dHev=540;                   % Latent heat of evaporation - [kcal/kg]
rhoL=1000;                  % Liquid Density - [kg/m^3]
rhoG=P/R/Tg0*MWair;         % Gas Density - [kg/m^3]
muG=2.3*1e-5;               % Gas Viscosity - [kg/m/s]
kG=8*1.0e-6;                % Gas Thermal Conductivity - [kcal/m/s/K]
Diff=1.8*1e-5;              % Diffusion Coefficient - [m^2/s]

% Antoine equation --- P = [mmHg] , T = [K]
A=18.3036;
B=3816.44;
C=-46.13;

% Process stream
Qmilk=(1750+10*D)/3600;             % Milk Flow Rate [kg/h -> kg/s]
fat=4.76/100;                % Fat Weight Fraction  [kg/kg]
Qfat=Qmilk*fat;              % Fat Flow Rate [kg/s]
Ql=Qmilk-Qfat;               % Liquid Flow Rate [kg/s]
Win=Ql/Qfat;                 % Inlet Moisture [kg/kg]
Wout=0.006;                  % Inlet Moisture [kg/kg]
dp=2*1e-4;                   % Particle Diameter [m]
vp0=0.25+0.01*D;                     % Initial particle velocity [m/s]
mp0=rhoL*3.14/6*(dp^3);      % Initial particle mass [kg]
mdry=mp0/(1+Win);            % Dry particle mass [kg]
mout=mdry*(1+Wout);          % Outlet particle mass [kg]
nAvp=Qmilk/mp0;              % [drops/s]

% Gas stream
Gdry=70e+3/3600;                    % Dry Gas Flow Rate [kg/s]
D=5 + 0.1*D;                      % Dryre Diameter [m]

% Momentum Balance
vg=Gdry/rhoG/3.14/D^2*4;                        % m/s
vs=(rhoL-rhoG)*dp^2*9.81/18/ muG;               % m/s
vp=vg+vs;                                       % m/s
vs0=vp0-vg;                                     % m/s

% Mass and heat transfer
Re=rhoG*vs*dp/muG;           % Reynolds
Pr=muG*CpG /kG;              % Prandtl
Sc= muG/rhoG/Diff;           % Schmidt
Nu=2+0.4*Re^0.5*Pr^(1/3);    % Nusselt
Sh=2+0.4*Re^0.5*Sc^(1/3);    % Sherwood
h=Nu*kG/dp;                  % kcal/m^2/s/K
Kc=Sh*Diff/dp;               % m/s
Kpw=Kc/R/Tg0*MWwater;        % kg/m^2/s/atm

% ODE Solver Setting
y0 = [mp0 0 Tp0 Tg0 vs0 0];
tspan = 0:0.1:10;   % [s]
[t,y] = ode15s(@dryer_model,tspan,y0);

Q5)

By increasing temperatures 

Opex: it decreases because we need smaller amount of steam in the heat exchanger

Capex: it alse decreases as the flux of air entering the dryer goes down. As a result, a smaller dryer is applied


% Graphics
figure(1)
plot(t,y(:,1)); xlabel('Time [s]'); ylabel('m_{P} [kg]'); grid on; hold off
figure(2)
plot(t,y(:,6)); xlabel('Time [s]'); ylabel('z [m]'); grid on 
figure(3)
plot(t,y(:,3)); hold on; plot(t,y(:,4)); xlabel('Time [s]'); ylabel('Temperature [K['); grid on 
legend('T_{P}','T_{G}','Location','best'); hold off
% Results
index = max((find(y(:,1)>2e-10)));
disp("at t = "+ t(index)+ ", z equals to "+y(index,6)+" m.")
function out = dryer_model(t,y)
global P Gdry A B C rhoL Kpw nAvp h MWwater MWair dHev CpL CpG rhoG muG dp mp0 vg

mp=y(1);
Gvap=y(2);
Tp=y(3);
Tg=y(4);
vs=y(5);
z=y(6);

Pw=MWair/MWwater*Gvap*P/Gdry;
P0=exp(A-B/(Tp+C))/760;                                                     % mmHg--->atm
Sp=3.14*((mp/rhoL*6/3.14)^(1/3))^2;

out(1)=Kpw*Sp*(Pw-P0);
out(2)=-Kpw*Sp*(Pw-P0)*nAvp;
out(3)=(h*Sp*(Tg-Tp)+Kpw*Sp*(Pw-P0)*dHev)/mp/CpL;
out(4)=-h*Sp*(Tg-Tp)*nAvp/(Gvap+Gdry)/CpG;
out(5)=((1-rhoG/rhoL)*9.81-3*vs*muG*3.14*dp/mp0);
out(6)=(vs+vg);

out=out';

end
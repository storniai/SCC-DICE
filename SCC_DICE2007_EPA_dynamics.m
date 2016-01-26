function [T_A,F,M_A,M_U,M_L] = SCC_DICE2007_EPA_dynamics(CS,ts,E,Fex,E_CH4,E_N2O);

% INPUTS:
%   CS  -- Climate sensitivity; eq temp change from CO2 doubling [deg C]
%   ts  -- Years defining time periods
%   E   -- Path of total CO2 emissions [GtC/decade]
%   Fex -- Path of non-CO2 forcings [W/m^2]
%   E_CH4 -- Path of CH4 emissions [MtC/decade]
%   E_N2O -- Path of N2O emissions [MtC/decade]
% OUTPUTS:
%   T_A -- Path of atmospheric temperature changes [deg C]
%   F   -- Path of total forcings [W/m^2]
%   M_A -- Path of atmospheric CO2 concentration [GtC]
%   M_U -- Path of upper ocean CO2 concentration [GtC]
%   M_L -- Path of lower ocean CO2 concentration [GtC]

years = ts;
M_A   = zeros(size(years)); % Atmospheric carbon concentration
M_Aav = zeros(size(years)); % 2-decade lag atmospheric carbon conc.
M_U   = zeros(size(years)); % Upper ocean carbon concentration
M_L   = zeros(size(years)); % Lower ocean carbon concentration
F     = zeros(size(years)); % Total forcing
T_A   = zeros(size(years)); % Atmospheric temperature
T_Aeq = zeros(size(years)); % Equilibrium atmospheric temp 
                                 % (conditional on current CO2 conc.)
T_L   = zeros(size(years)); % Lower ocean temperature
CH4ppb= zeros(size(years));
N2Oppv= zeros(size(years));

% DICE2007 PARAMETERS:
% Carbon transport coefficients:
b = [0.810712 0.189288 0.000000 ; ...
     0.097213 0.852787 0.050000 ; ...
     0.000000 0.003119 0.996881];
%b(1,1) = 1 - b(1,2);
%b(2,1) = 587.473*b(1,2)/1143.894;
%b(2,2) = 1 - b(2,1) - b(2,3);
%b(3,2) = 1143.894*b(2,3)/18340;
%b(3,3) = 1 - b(3,2);
% Estimated forcings of equilibrium CO2 doubling: 
FCO22X = 3.8;
% Atmospheric temperature dynamics coefficients:
C = [0.220; 0; 0.3; 0.05];
% Starting values for concentrations in 2005 [GtC]:
M_A(1) = 808.9;  
M_U(1) = 1255;   
M_L(1) = 18365;  


% Carbon concentration transition equations:
t = 1;
for year = years(2:end)'; t = t + 1;
    M_A(t) = b(1,1)*M_A(t-1) + b(2,1)*M_U(t-1) + E(t-1);
    M_U(t) = b(1,2)*M_A(t-1) + b(2,2)*M_U(t-1) + b(3,2)*M_L(t-1);
    M_L(t) = b(2,3)*M_U(t-1) + b(3,3)*M_L(t-1);
end;
t = 0;
for year = years(1:end-1)'; t = t + 1;
    M_Aav(t) = (M_A(t)+M_A(t+1))/2;
end;
M_Aav(end) = 0.9796*M_Aav(end-1);


% DICE2007 default exogenous non-CO2 forcings:
for zz = 1:1; break; 
    Fex  = zeros(size(years)); 
    Fex0 = -.06;     % Estimate of 2000 forcings of non-CO2 GHGs
    Fex1 = 0.3;      % Estimate of 2100 forcings of non-CO2 GHGs 
    for t = 1:length(years);
        Fex(t) = Fex0 + 0.1*(Fex1-Fex0)*(t-1)*(t<12) + 0.36*(t>=12);
    end;
end;


T_A(1)  = 0.7307; % Atmospheric temp change from 1900 to 2000 [deg C]
T_Aeq(1)= 0.7307;
T_L(1)  = 0.0068; % Lower ocean strata temp change 1900 to 2000 [deg C]
t = 1;
for year = years(2:end)'; t = t + 1;
    F(t)      = FCO22X * log((M_Aav(t)+0.000001)/596.4)/log(2) + Fex(t);   
    T_A(t)    = T_A(t-1) + C(1)*( F(t) - FCO22X/CS*T_A(t-1) - ...
                C(3)*( T_A(t-1) - T_L(t-1) ) );
    T_Aeq(t)  = F(t)*CS/FCO22X;
    % The following line was added to avoid over- and under-shoots to eq T
    if CS < .5; T_A(t) = T_Aeq(t); end; 
    T_L(t) = T_L(t-1) + C(4)*( T_A(t-1) - T_L(t-1) );
end;
   
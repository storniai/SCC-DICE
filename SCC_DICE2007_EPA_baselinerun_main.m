%==========================================================================
% SCC_DICE2007_EPA_baselinerun_damageTFP.m
%
% This file modifies the original EPA SCC code to investigate the
% sensitivity of SCC to alternative damages specifications
%
% The modification to the original code are:
%   uses multiple discount rates 1%...10%
%   uses climate sensitivity of Roe & Baker distribution function
%       implemented by Kopp, R.
%   modified output printing (.xls)
%==========================================================================

% ADDED: parameter specifying share of damages to TFP
phi = 0;
% ADDED: define alternative damage parameters
a3 = ((1-.15)^(-1)-1)/2.5^2; %calibrated to be 15% 
a4 = ((1-.30)^(-1)-1)/2.5^2; %calibrated to be 30% 

% Input parameters:
rate = [.01 .02 .03 .04 .05 .06 .07 .08 .09 .10];    % Fixed discount rates
a    = 0.0028388; % Coefficient of damage function
b    = 2;         % Exponent of damage function
H    = 2300;      % Time horizon for calculating the SCC [year]

% Climate sensitivity parameters TO BE UPDATED MORE
fbar = 0.61979;
fsig = 0.18407;
maxT = 10;
cs0 = 1.2;

%Monte Carlo simulation
MC = 1000;
iterations = 5*3*MC; count = 0; starttime = cputime;

for scenario = 1:5;
results = zeros(length(rate),41); % 41 is nrof 10-year steps
for index = 1:length(rate);        
%--------------------------------------------------------------------------
% Import data:
% NOTES: GDP input data should be in trillions of US$2005
%        Population input data should be in millions of individuals
%        CO2 emissions input data should be in GtC/decade
%        nonCO2 forcings input data should be in W/m^2
InputFile = [pwd '/SCC_DICE2007_EPA_input.xls'];
Y    = xlsread(InputFile,'GDP'); EMFyears = Y(:,1); 
N    = xlsread(InputFile,'Population');      
E    = xlsread(InputFile,'IndustrialCO2');  
El   = xlsread(InputFile,'LandCO2');    
Fex1 = xlsread(InputFile,'EMFnonCO2forcings'); 
Fex2 = xlsread(InputFile,'OthernonCO2forcings'); 
inflator = 122.58/114.52; % World GDP inflator 2007/2005
if scenario == 1;     scenarioname = 'IMAGE';
    Y=Y(:,2)*inflator;N=N(:,2);E=E(:,2)+El(:,2);Fex=Fex1(:,2)+Fex2(:,5);
elseif scenario == 2; scenarioname = 'MERGEoptimistic'; 
    Y=Y(:,3)*inflator;N=N(:,3);E=E(:,3)+El(:,3);Fex=Fex1(:,3)+Fex2(:,5); 
elseif scenario == 3; scenarioname = 'MESSAGE';
    Y=Y(:,4)*inflator;N=N(:,4);E=E(:,4)+El(:,4);Fex=Fex1(:,4)+Fex2(:,5); 
elseif scenario == 4; scenarioname = 'MiniCAMbase';
    Y=Y(:,5)*inflator;N=N(:,5);E=E(:,5)+El(:,5);Fex=Fex1(:,5)+Fex2(:,5);
elseif scenario == 5; scenarioname = '5thScenario';
    Y=Y(:,6)*inflator;N=N(:,6);E=E(:,6)+El(:,6);Fex=Fex1(:,6)+Fex2(:,5); 
end;
% Interpolate EMF paths to DICE2007 time periods
years = [2005:10:2405]';
warning('off','all');
Y = interp1(EMFyears,Y,years(1:length(EMFyears)),'linear')'; Y(isnan(Y))=0;
temp = zeros(size(years')); temp(1:length(Y)) = Y; Y = temp';
N = interp1(EMFyears,N,years(1:length(EMFyears)),'linear')'; N(isnan(N))=0;
temp = zeros(size(years')); temp(1:length(N)) = N; N = temp';
E = interp1(EMFyears,E,years(1:length(EMFyears)),'linear')'; E(isnan(E))=0;
temp = zeros(size(years')); temp(1:length(E)) = E; E = temp';
Fex = interp1(EMFyears,Fex,years(1:length(EMFyears)),'linear')'; Fex(isnan(Fex))=0;
temp = zeros(size(years')); temp(1:length(Fex)) = Fex; Fex = temp';
warning('on','all');
y = Y./(N+eps);
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Solve for implied path of exogenous technical change
A     = zeros(length(years),1);
A(1)  = 0.02722; % From DICE2007
gamma = 0.3;     % From DICE2007
delta = 0.1;     % From DICE2007
s     = 0.22;    % Approximate optimal savings in DICE2007, here fixed
K(1)  = (Y(1)/A(1)/(N(1)^(1-gamma)))^(1/gamma);
for t = 2:length(years);
    K(t) = K(t-1)*(1-delta)^10 + s*Y(t-1)*10;
    A(t) = Y(t)/(N(t)+eps)^(1-gamma)/(K(t)+eps)^gamma;
end;

% ADDED: solve for exogenous TFP growth rate
Arate = zeros(length(years),1);
for t = 1:length(years)-1;
    Arate(t) = A(t+1)/A(t) - 1;
end;
%--------------------------------------------------------------------------

%Monte Carlo simulation
SCC_MC = zeros(MC,length(years)); 
for mc = 1:MC;
    
    % Take random draw of CS:
    CS = icdfRoeBaker(rand,fbar,fsig,maxT,cs0)
               
    % Calculate SCC(t) for each time period
    L0  = zeros(length(years),1);
    L1  = zeros(length(years),1);
    SCC = zeros(length(years),1);

    % Reference temp path and losses:
    E0 = E;
    [T0,F,M_A,M_U,M_L] = SCC_DICE2007_EPA_dynamics(CS,years,E0,Fex); 
    L0 = 1-1./(1+a*T0.^b); 
    L0A = phi.*L0;            %share of damages to TFP (A)
    
    if phi==0; % if no damages to TFP
         L0Y = L0;
    else L0Y = 1- (1-L0)./(1-L0A); %share on output; defined so that product of (1-L0A)*(1-L0Y)=(1-L0)
    end;
    
    % ADDED:
     A0 = zeros(length(years),1); %initialize endogenous TFP arrays
     A1 = zeros(length(years),1); 
        A0(1)=A(1);
        A1(1)=A(1);
    
    % Reference GDP and consumption paths:
    Y0 = zeros(length(years),1); 
    K0 = zeros(length(years),1); K0(1) = K(1);
    C0 = zeros(length(years),1);
    for tau = 1:length(years(find(years<=EMFyears(end))));
        Y0(tau)   = A0(tau)*K0(tau)^gamma*N(tau)^(1-gamma)*(1-L0Y(tau));
        C0(tau)   = (1-s)*Y0(tau);
        K0(tau+1) = K0(tau)*(1-delta)^10 + s*10*Y0(tau);
        A0(tau+1) = A0(tau)*(1+Arate(tau))*(1-L0A(tau))^10;
    end;     
                    
    % SCC as change in consumption discounted at a fixed rate:
    t = 0;
    for year = years(find(years<=EMFyears(end)))'; t = t + 1;        
             
        % Perturbed temp path and losses:
        E1 = E; E1(t) = E1(t)+1;
        [T1,F,M_A,M_U,M_L] = SCC_DICE2007_EPA_dynamics(CS,years,E1,Fex);
        L1 = 1-1./(1+a*T1.^b);
        L1A = phi.*L1;            %share of damages to TFP (A)
    
    if phi==0; % if no damages to TFP
         L1Y = L1;
    else L1Y = 1- (1-L1)./(1-L1A); %share on output; defined so that product of (1-L0A)*(1-L0Y)=(1-L0)
    end;
        
        % Perturbed GDP and consumption paths:
        Y1 = Y0; K1 = K0; C1 = C0;
        for tau = t:length(years(find(years<=EMFyears(end))));
            Y1(tau)   = A1(tau)*K1(tau)^gamma*N(tau)^(1-gamma)*(1-L1Y(tau));
            C1(tau)   = (1-s)*Y1(tau);
            K1(tau+1) = K1(tau)*(1-delta)^10 + s*Y1(tau)*10;
            A1(tau+1) = A1(tau)*(1+Arate(tau))*(1-L1A(tau))^10;
        end;
            
        % Calculate SCC: 
        C0i = interp1(years(t:end),C0(t:end),years(t):years(end),'linear')'; 
        Ni  = interp1(years(t:end),N(t:end),years(t):years(end),'linear')'; 
        C1i = interp1(years(t:end),C1(t:end),years(t):years(end),'linear')'; 
        T0i = interp1(years(t:end),T0(t:end),years(t):years(end),'linear')'; 
                
        % Constant discount rate:
        for zz = 1:1; %break;
        tt = 0;
        for tau = years(t):H-1; tt = tt + 1;
            SCC(t) = SCC(t) + ( C0i(tt) - C1i(tt) ) * 1/((1+rate(index))^(tt-1)); 
        end;
        end;
        
	end;
    
    % ADDED - computation of output per capita
        y  = 10^6 * Y ./ (N+eps);  %GDP is in 10^12; pop in 10^6
        y0 = 10^6 * Y0 ./ (N+eps); %GDP is in 10^12; pop in 10^6 
    %calculate scaled GDP as share of t=1 GDP, and scaled TFP 
    %initiate vectors
        y_p  = zeros(length(years),1);
        y0_p = zeros(length(years),1);
        A_p  = zeros(length(years),1);
    %loop through years
        for i=1:length(years)
            y_p(i)  = y(i)  / y(1);
            y0_p(i) = y0(i) / y0(1);
            A_p(i) = A(i) / A(1);
        end;
    
    % Convert SCC path from trillion$/GtC/yr to $/tCO2/yr:
    SCC = SCC / (10^-3) * 12/44; 
    SCC_MC(mc,:) = SCC'; %building Monte Carlo summary
end; %Monte Carlo run

results(index,:) = mean(SCC_MC);
end;  %rate loop

header = [{'year/rate'},cellstr(num2str(rate(:)))'];
data = [years,results'];
sheet = [scenarioname];
xlswrite([pwd ['\SCC_DICE2007_EPA_output_avg_' date '.xls']],[header;num2cell(data)],sheet);
end;  %scenario loop
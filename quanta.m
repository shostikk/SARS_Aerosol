function [n,R,N,N_in,ER,R0] = quanta(time,time_I,V,I,AER,dt_in,IRs,EM,TR)
%QUANTA Infection model 
%   Model expressed Buonanno et al. (2020) Environment International 141,105794
%   https://doi.org/10.1016/j.envint.2020.105794
%
%   [n,R,N,Nin,ER,R0] = quanta(time,time_I,V,I,AER,IR,EM,TR)
%
%   n          number of quanta (quanta m^(-3)): a quantum is defined
%	       as the dose of airborne droplet nuclei required to cause infection in 63%
%              of susceptible persons
%   R          individual risk
%   N          number of droplets per m^3
%   N_in       number of inhaled droplets for a target person arriving at 
%              all values of time
%   ER         quanta emission rate (Eq 2)
%   R0         Basic reproduction number
%
%   time       time vector (h)
%   time_I     duration of visit for infected persons (h), starts from zero
%   V          slace volume (m3)
%   I          number of infected persons
%   AER        air exhance rate of ventilation (h^-1)
%   dt_in      duration of visit for target persons (h), calculated for all
%              time -vector values
%   IR         Inhalation rate (0=rest, 1=stand, 2=light, 3=moderate,
%                               4=heavy, 5=mean(stand,light)
%   EM         Emission model (0=voiced, 1=whispered, 2=vocal,
%                              3=breath, 4=speak, 5=average)
%   TR         (optional) rate of target people entrance h^-1)
%
%   The quanta density calcuation is integrated numerically, instead of
%   using Eq 3.
%
% Coded June 1, 2020, by Simo Hostikka

if (nargin < 9), TR = 0;end
%
dt = time(2)-time(1);
Ntot = length(time);
%
% Parameters
c_i = 0.02; 	% Ratio between one infectious quantum and the infectious dose (viral NRA copies)
c_v = 10^9; 	% RNA copies per mL
c_v = c_v*10^6; % per m3 of droplet
% Droplet size groups (micron)
D = [0.8 1.8 3.5 5.5];
V_i = pi*(4/3)*((D*10^(-6))/2).^3; % droplet volumes m3
%
% Particle concentrations for size groups, part/cm3
N_i_voiced = [0.236 0.068 0.007 0.011];
N_i_whispered = [0.11 0.014 0.004 0.002];
N_i_vocal = [0.751 0.139 0.139 0.059];
N_i_breath = [0.084 0.009 0.003 0.002];
N_i_speak = 0.5*(N_i_vocal+N_i_voiced);
N_i_ave = 0.25*(N_i_voiced+N_i_whispered+N_i_vocal+N_i_breath);
%
switch EM
    case 0
        N_i = N_i_voiced;
    case 1
        N_i = N_i_whispered;
    case 2
        N_i = N_i_vocal;
    case 3
        N_i = N_i_breath;
    case 4
        N_i = N_i_speak;
    case 5
        N_i = N_i_ave;
end
N_i = N_i*(100^3); % cm3 to m3
%
% Exhalation rate m3/h
IR_rest = 0.49;
IR_stand = 0.54;
IR_light =  1.38;
IR_moderate = 2.35;
IR_heavy = 3.30;
%
% CHOOSE HERE
switch IRs
    case 0
        IR = IR_rest;
    case 1
        IR = IR_stand;
    case 2
        IR = IR_light;
    case 3
        IR = IR_moderate;
    case 4
        IR = IR_heavy;
    case 5
        IR = 0.5*(IR_stand+IR_light);
end
IR_target = IR;
%
%
% Calculations
%
ER = c_i*c_v*IR*sum(N_i.*V_i) %quanta / h
NR = IR*sum(N_i); % droplets/h
%
%
n0 = 0;
k = 0.24; % h^-1
lambda = 0.63; % h^-1
IVRR = AER + k + lambda;
n(1) = n0;
N(1) = 0;
Nsteps_in = round(dt_in/dt);
time2 = min(time):dt:(Ntot+Nsteps_in)*dt;
% calculate number of droplets and quanta per m3
for i = 2:(Ntot+Nsteps_in)
    if (and((time2(i)<= time_I),(time2(i)>=0)))
        ERq = ER;
        NRq = NR;
    else
        ERq = 0;
        NRq = 0;
    end
    n(i) = n(i-1) + dt*((ERq*I/V) - n(i-1)*IVRR);
    N(i) = N(i-1) + dt*((NRq*I/V) - N(i-1)*IVRR);
end
%
% Dose integration loop
R0 = 0;
for i = 1:Ntot
    n_in = IR_target*dt*sum(n(i:i+Nsteps_in));
    R(i)=1-exp(-n_in);
    R0 = R0 + R(i)*dt*TR;
    Ntmp = IR_target*dt*sum(N(i:i+Nsteps_in));
    N_in(i)=Ntmp;
end
n = n(1:Ntot);
N = N(1:Ntot);
N_in = N_in(1:Ntot);
end

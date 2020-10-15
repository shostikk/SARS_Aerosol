% Example file for using quanta.m to calculate the concentrations and infection risk
% in a public bus.
% With small changes, this file can be used to study other spaces where the assumptions hold.

% time inside the bus (h)
dt_in = 10/60;
% time vector (h)
t = -dt_in:1/3600:2;
% bus volume (m3)
V = 75;
% number of passengers entering the bus per hour
TR = 180;
% air exchange rate (h^-1)
AER = 1;
% Switch for exhalaling/inhalation rate 
% (0=rest, 1=stand, 2=light, 3=moderate, 4=heavy, 5=mean(stand,light)
IR = 1;
% Switch for emission model
% (0=voiced, 1=whispered, 2=vocal,3=breath, 4=speak, 5=average)
EM=4;
% Base case is 10 min travel, speaking person, 1 h-1 exchange rate
[n_1,R_1,N_1,Nintmp,ER,R0_1] = quanta(t,dt_in,V,1,AER,dt_in,IR,EM,TR); 
%
%
% Effect of the air exchange rate
AER = 10;
[n_10,R_10,N_10,Nintmp,ER,R0_10] = quanta(t,dt_in,V,1,AER,dt_in,IR,EM,TR); 
figure(1)
[ax,h1,h2]=plotyy(t*60,[N_1; N_10],t*60,[R_1; R_10]*100);
R_max = 1; dR = 0.25;
Rstr='R_{10 min}(%)';
figfeat(ax,h1,h2,dt_in,R_max,dR,Rstr)
legend('N (1 h^{-1})','N (10 h^{-1})',...
    ['R (1 h^{-1}) R_0=' num2str(R0_1,2)],...
    ['R (10 h^{-1}) R_0=' num2str(R0_10,2)])
title('V = 75 m^3, 180 pax/h, standing, speaking')
figsaver('Bus_AER')
%
%
% Effect of the travelling time
AER = 1;
[n_1_5,R_1_5,N_1_5,Nintmp,ER,R0_1_5] = quanta(t,dt_in*0.5,V,1,AER,dt_in*0.5,IR,EM,TR);
figure(2)
[ax,h1,h2]=plotyy(t*60,[N_1; N_1_5],t*60,[R_1; R_1_5]*100);
R_max = 1; dR = 0.25;
Rstr='R_{10 min, 5 min}(%)';
figfeat(ax,h1,h2,dt_in,R_max,dR,Rstr)
legend('N (10 min)','N (5 min)',...
    ['R (10 min) R_0=' num2str(R0_1,2)],...
    ['R (5 min) R_0=' num2str(R0_1_5,2)])
title('V = 75 m^3, 180 pax/h, standing, speaking, 1 h^{-1}')
figsaver('Bus_dt_in')
%
%
% Effect of emission model
figure(3)
AER = 1;
EM = 3; % breath
[n_1b,R_1b,N_1b,Nintmp,ER,R0_1b] = quanta(t,dt_in,V,1,AER,dt_in,IR,EM,TR); 
[ax,h1,h2]=plotyy(t*60,[N_1; N_1b],t*60,[R_1; R_1b]*100);
R_max = 1; dR = 0.25;
Rstr='R_{10 min} (%)';
figfeat(ax,h1,h2,dt_in,R_max,dR,Rstr)
legend('N (speak)','N (breath)',...
    ['R (speak) R_0=' num2str(R0_1,2)],...
    ['R (breath) R_0=' num2str(R0_1b,2)])
title('V = 75 m^3, 180 pax/h, standing, 1 h^{-1}')
figsaver('Bus_EM_mode')

%
function figfeat(ax,h1,h2,dt_in,R_max,dR,Rstr)
set(ax(1),'XLim',60*[-dt_in 1])
set(ax(2),'XLim',60*[-dt_in 1])
set(ax(2),'YLim',[0 R_max],'YTick',[0:dR:R_max])
set(ax(1),'FontSize',14)
set(ax(2),'FontSize',14)
set(h1,'LineWidth',1.6)
set(h2,'LineWidth',1.6,'LineStyle','--')
xlabel('Entrance time (min)')
set(ax(1).YLabel,'String','N (part. m^{-3})')
set(ax(2).YLabel,'String',Rstr);
set(ax(1),'Position', [0.1458    0.1452    0.7083    0.7798])
set(ax(2),'Position', [0.1458    0.1452    0.7083    0.7798])

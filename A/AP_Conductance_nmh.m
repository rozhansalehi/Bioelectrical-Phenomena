 %% Hodgkin-Huxley model using Euler's method

%This script allows up to two current pulses to be injected into the membrane to excite a response. The widths can be independently set, but the pulses cannot overlap. If pulse 2 is turned on before pulse 1 ends it replaces pulse 1.

%% Parameteres
clc
clear 
close all
%maximal conductances (mS/cm^2); 1=K, 2=Na, 3=leakage;
 g(1)=36; g(2)=120; g(3)=0.3;
 
% membrane capacitance is 1 uF/cm^2
 Cm=1;
 
% battery voltage (mV) relative to resting potential; 1=K, 2=Na; 3=lk
 E(1)=-12; E(2)=115; E(3)=10.613;
 

%% variable initialization
%I_ext is in microamps/cm^2
 I_ext=0; V=-10; x=zeros(1,3); x(3)=1; t_rec=0;
 

%% applied pulses parameters
t_final=120; % t_final sets the time span of the simulation
I_on1=2.5; % I_on1 is the amplitude of pulse 1 in microamps/square cm
I_on2=2.5; % I_on2 is the amplitude of pulse 2 in microamps/square cm
T1=10; T2=28;  %T1=on time for pulse 1, T2=on time for pulse 2
Tw1=5; Tw2=5;  % Tw1 is the width of pulse 1, Tw2 is the width of pulse 2
dt=0.01;  % time step for integration in milliseconds
  
% All times are in milliseconds. If T2>t_final the second pulse is not applied.

%% integration by Euler's method

% computations for t < 0 establish initial conditions at t = 0.

         for t=-30:dt:t_final
             if t==T1; I_ext=I_on1; end %turns on external current at t=T1
             if t==T1+Tw1; I_ext=0; end %turns off external current
             if t==T2; I_ext=I_on2; end %turns on second pulse
             if t==T2+Tw2; I_ext=0; end %turns off second pulse
     
         %alpha parameters in H-H model
         alpha(1)=(10-V)/(100*(exp((10-V)/10)-1));
         alpha(2)=(25-V)/(10*(exp((25-V)/10)-1));
         alpha(3)=0.07*exp(-V/20);
 
        %beta parameters in H-H model
         beta(1)=0.125*exp(-V/80);
         beta(2)=4*exp(-V/18);
         beta(3)=1/(exp((30-V)/10)+1);
 
         %time constants (msec) and asymptotic values
         tau=1./(alpha+beta);
         x_0=alpha.*tau;
    
        %Euler integration
        x=(1-dt./tau).*x + dt./tau.*x_0;
 
         %conductance calculations
        gnmh(1)=g(1)*x(1)^4;
        gnmh(2)=g(2)*x(2)^3*x(3);
        gnmh(3)=g(3);
 
       %membrane voltage update
        I=gnmh.*(V-E);
        V=V+(1/Cm)*dt*(I_ext-sum(I));
 
        %plotting records
       
       
        if t>=0
            t_rec=t_rec+1;
            x_axis(t_rec)=t;
            y_axis(t_rec)=V;
            G(t_rec,1)=gnmh(1); % GK
            G(t_rec,2)=gnmh(2); % GNa
            n(t_rec)=x(1);
            m(t_rec)=x(2);
            h(t_rec)=x(3);
        end
           end
        
 
plot(x_axis,y_axis,'c','linewidth',1.5); 
title('Action Potential')
xlabel('Time (ms)'); 
ylabel('Relative Membrane Voltage (mV)');
 

 figure;
 plot(x_axis,G,'linewidth',1.5);
 title('Ion Conductance')
 xlabel('Time (ms)');
 ylabel('Conductance (mS/cm^2)');
 legend('K','Na')
 
 figure;hold on;
 plot(x_axis,n,'color',[1 0 0.4],'Linewidth',1.5); 
 plot(x_axis,m,'color',[120 50 150]/255,'Linewidth',1.5);
 plot(x_axis,h,'color',[0 1 0.1],'Linewidth',1.5); xlabel('Time (ms)'); ylabel('n,m,h(t)');
 legend('n', 'm','h')


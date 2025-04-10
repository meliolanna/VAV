
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   VIBRATION ANALYSIS AND VIBROACOUSTICS   %
%   %   %   %   %   %   %   %   %   %   %   %
%               ASSIGNMENT 1                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% DATA
% body 1
m1 = 5;     % mass, [kg]
J1 = 10;    % moment of inertia, [kgm^2]
Ra = 0.80;  % radius of the largest disk, [m]
Rb = 0.50;  % radius of the ineer disk, [m]
% body 2
m2 = 7;     % mass, [kg]
J2 = 8;     % moment of inertia, [kgm^2]
R2 = 0.60;  % radius, [m]
% spring and dampers
k1 = 2000;  % elastic coeff, [N/m]
c1 = 10;    % damping coef, [Ns/m]
k2 = 2300;  % elastic coeff, [N/m]
c2 = 35;    % damping coef, [Ns/m]
k3 = 2500;  % elastic coeff, [N/m]
c3 = 20;    % damping coef, [Ns/m]
%parameters
A0 = 200;   % [Nm]
f0 = 1;     % [Hz]
phi0 = pi/3;
Ai = 200;
fi = [0.5, 20];
B1 = 400;
f1 = 0.65;
phi1 = -pi/2;
B2 = 44.4;
f2 = 1.95;
phi2 = pi/2;
B3 = 16;
f3 = 3.25;
phi3 = -pi/2;

% 1.a
Jeq = (m1*(R2^2))/4 + J1*(R2/(2*Ra))^2 + J2;
keq = k1*(Ra + Rb)^2*(R2/(2*Ra))^2 + k2*(R2^2)/4 + k3*(R2^2);
ceq = c1*(Ra + Rb)^2 * (R2/(2*Ra))^2 + c2*(R2^2)/4 + c3*(R2^2);
% 1.b
lambda = roots([Jeq ceq keq]); % solutions of the characteristic eq
omega_0 = sqrt(keq/Jeq); % natural frequency
% 1.c
ccr = 2*Jeq*omega_0; % critical damping
h = ceq/ccr; % non dimentional damping ratio
omega_d = omega_0 * sqrt(1-h^2); % damped frequency
alpha = h*omega_0;

% FREE MOTION
% 2.a
theta_20 = 0.1; % IC arbitrary chosen values for initial conditions
omega_20 = 5;   % IC same
t = 0:0.01:8; % time [s]

A = theta_20;
B = (omega_20 + alpha*theta_20)/omega_d;
theta_2 = exp(-alpha*t) .* (A*cos(omega_d*t) + B*sin(omega_d*t));

figure(1)
plot(t,theta_2,'r','LineWidth',1.2); 
grid minor
xlabel('Time [s]','Fontsize',10); 
ylabel('Angular displacement [rad]','Fontsize',10);
%legend('Sinusoidal','Exponential','Fontsize',9)
title('Time response','Fontsize',12);

% 2.b (using 5h)
hb = 5*h;
omega_db = omega_0 * sqrt(1-hb^2); 
alphab = hb*omega_0;
theta_2b = exp(-alphab*t) .* (A*cos(omega_db*t) + B*sin(omega_db*t));

figure(2)
plot(t,theta_2b,'b','LineWidth',1.2); 
grid minor
xlabel('Time [s]','Fontsize',10); 
ylabel('Angular displacement [rad]','Fontsize',10);
%legend('Sinusoidal','Exponential','Fontsize',9)
title('Time response using 5h','Fontsize',12);

% 2.c (using 25h) SBaaaaagliato
hc = 25*h;
%omega_dc = omega_0 * sqrt(1-hc^2); 
alphac = hc*omega_0;
L1 = -alphac + sqrt(alphac^2 - omega_0^2);
L2 = -alphac - sqrt(alphac^2 - omega_0^2);
T1 = (theta_20*L2 - omega_20)/(L2-L1);
T2 = (theta_20*L1 - omega_20)/(L1 - L2);
%theta_2c = exp(-alphac*t) .* (A*cos(omega_dc*t) + B*sin(omega_dc*t));
theta_2c = T1*exp(L1*t) + T2*exp(L2*t);

figure(3)
plot(t,theta_2c,'g','LineWidth',1.2); 
grid minor
xlabel('Time [s]','Fontsize',10); 
ylabel('Angular displacement [rad]','Fontsize',10);
%legend('Sinusoidal','Exponential','Fontsize',9)
title('Time response using 25h','Fontsize',12);

figure(7)
plot(t, theta_2, 'r')
grid minor
hold on
plot(t, theta_2b, 'b')
hold on
plot(t, theta_2c, 'g')
xlabel('Time [s]', 'FontSize',10)
ylabel('Angular displacement [rad]', 'FontSize', 10)

% FORCED MOTION
% 3.a 
freq = linspace(0,50, 500); %frequenze per rappresentare
% h
FRF = 1./(-Jeq*freq.^2 - ceq*freq.*1i + keq);
mod_FRF = abs(FRF);
phase_FRF = angle(FRF);
% 5h
FRF_5 = (omega_0^2 - freq.^2 + 2*1i*omega_0*hb.*freq)./(Jeq*(omega_0^2-freq.^2).^2 + (2*omega_0*hb.*freq).^2);
mod_FRF5 = abs(FRF_5);
phase_FRF5 = angle(FRF_5);
% 25h
FRF_25 = (omega_0^2 - freq.^2 + 2*1i*omega_0*hc.*freq)./(Jeq*(omega_0^2-freq.^2).^2 + (2*omega_0*hc.*freq).^2);
mod_FRF25 = abs(FRF_25);
phase_FRF25 = angle(FRF_25);

%questi grafici fanno tremendamente schifo ma intanto li abbozzo
figure(4)
subplot(321)
plot(freq, mod_FRF);
grid minor
xlabel("\Omega [rad/s]");
ylabel("|H(\Omega)| [rad/N]");
subplot(322)
plot(freq, phase_FRF);
grid minor
xlabel("\Omega [rad/s]");
ylabel("Phase H(\Omega) rad")

subplot(323)
plot(freq, mod_FRF5);
grid minor
xlabel("\Omega [rad/s]");
ylabel("|H5(\Omega)| [rad/N]");
subplot(324)
plot(freq, phase_FRF5);
grid minor
xlabel("\Omega [rad/s]");
ylabel("Phase H5(\Omega) rad")

subplot(325) 
plot(freq, mod_FRF25);
grid minor
xlabel("\Omega [rad/s]");
ylabel("|H25(\Omega)| [rad/N]");
subplot(326)
plot(freq, phase_FRF25);
grid minor
xlabel("\Omega [rad/s]");
ylabel("Phase H25(\Omega) rad")


% 3.b
Omega = 2*pi*1; %valore a caso
phi = pi/3; % valore a caso

FRF_Omega = 1/(-Jeq*Omega^2 - ceq*Omega*1i + keq);
theta_p = abs(FRF_Omega) * A * cos(Omega*t + phi + angle(FRF_Omega)*1i); 
theta_tr = theta_2 + theta_p;

figure(5)
plot(t,theta_tr); 
grid minor
xlabel('Time [s]'); 
ylabel('Angular displacement [rad]');
title('Time response to the torque','Fontsize',12);

% 3.c
Omegai1 = 2*pi*0.5;
Omegai2 = 2*pi*20;

FRF_Oi1 = 1/(-Jeq*Omegai1^2 - ceq*Omegai1*1j + keq);
theta_tri1 = theta_2 + abs(FRF_Oi1) * A * cos(Omegai1*t + angle(FRF_Oi1)*1j); 
FRF_Oi2 = 1/(-Jeq*Omegai2^2 - ceq*Omegai2*1j + keq);
theta_tri2 = theta_2 + abs(FRF_Oi2) * A * cos(Omegai2*t + angle(FRF_Oi2)*1j); 

figure(6) %Spoiler: vengono ancora tutti uguali ma non dovrebbero
subplot(1,2,1)
plot(t, theta_tri1)
grid minor
xlabel('Time [s]'); 
ylabel('Angular displacement [rad]');
subplot(1,2,2)
plot(t, theta_tri2)
grid minor
xlabel('Time [s]'); 
ylabel('Angular displacement [rad]');

%3.d

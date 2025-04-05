
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

% 2.c (using 25h)
hc = 25*h;
omega_dc = omega_0 * sqrt(1-hc^2); 
alphac = hc*omega_0;
theta_2c = exp(-alphac*t) .* (A*cos(omega_dc*t) + B*sin(omega_dc*t));

figure(3)
plot(t,theta_2c,'g','LineWidth',1.2); 
grid minor
xlabel('Time [s]','Fontsize',10); 
ylabel('Angular displacement [rad]','Fontsize',10);
%legend('Sinusoidal','Exponential','Fontsize',9)
title('Time response using 25h','Fontsize',12);

% FORCED MOTION
% 3.a 
Omega = 1; %valore a caso
phi = 3; % valore a caso

% h
FRF = 1/(-Jeq*Omega^2 - ceq*Omega*j + keq);
mod_FRF = abs(FRF);
phase_FRF = angle(FRF);
% 5h
FRF_5 = (omega_0^2 - Omega^2 + 2*j*omega_0*hb*Omega)/(Jeq*(omega_0^2-Omega^2)^2 + (2*omega_0*hb*Omega)^2);
mod_FRF5 = abs(FRF_5);
phase_FRF5 = angle(FRF_5);
% 25h
FRF_25 = (omega_0^2 - Omega^2 + 2*j*omega_0*hc*Omega)/(Jeq*(omega_0^2-Omega^2)^2 + (2*omega_0*hc*Omega)^2);
mod_FRF25 = abs(FRF_25);
phase_FRF25 = angle(FRF_25);

%theta_p = mod_FRF * C_0 * cos(Omega*t + phi + angle(FRF)*j); 

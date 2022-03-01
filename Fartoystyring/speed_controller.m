%% Information 
% This file is only an example of how you can start the simulation. The
% sampling time decides how often you store states. The execution  time
% will increase if you reduce the sampling time.

% Please note that the file "pathplotter.m" (only used in the second part
% of the assignment) shows the ship path during the path following and
% target tracking part of the assignment. It can be clever to adjust the sampling
% time when you use that file because it draws a sketch of the ship in the
% North-East plane at each time instant. Having a small sampling time will
% lead to multiple ship drawings on top of each other. 

% You should base all of your simulink models on the MSFartoystyring model
% and extend that as you solve the assignment. For your own sake, it is
% wise to create a new model and run file for each task.

% The msfartoystyring.m file includes the ship model. You are not allowed
% to change anything within that file. You need to include that file in
% every folder where you have a simulink model based on
% "MSFartoystyring.slx". 

% WP.mat is a set of six waypoints that you need to use in the second part of
% the assignment. The north position is given in the first row and the east
% position in the second row. 

%% System information
L_pp = 304.8; % [m]
delta_max = deg2rad(25); % [deg]
n_max = (85*2*pi)/60; % [rad/s]

%% Simulation
tstart=0;           % Sim start time
tstop=10000;        % Sim stop time
tsamp=1;           % Sampling time for how often states are stored. (NOT ODE solver time step)
                
p0=zeros(2,1);      % Initial position (NED)
v0=[3 0]';       % Initial velocity (body)
psi0=0;             % Inital yaw angle
r0=0;               % Inital yaw rate
c=1;                % Current on (1)/off (0)

sim speed_control % The measurements from the simulink model are automatically written to the workspace.


figure(1)
title('Control Performance of Velocity')
plot(t, u_d(1: 10001), 'b');
hold on
plot(t, v(:, 1), 'r');
axis([0 10000 0 10]);
legend('Referece u', 'Actual u')
xlabel('Time (s)')
ylabel('u (m/s)')

grid on

figure(2)
title('Control Performance of Psi')
plot(t, psi_d, 'b');
plot(t, psi, 'r');
hold on

xlabel('Time (s)')
ylabel('Psi (rad)')
axis([0 10000 -0.003 0.001]);
line([0 10000], [0 0]);
legend( 'Actual Psi', 'Referece Psi')
grid on













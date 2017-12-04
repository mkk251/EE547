% Kyle Phan and Mulugeta Kembata (1368546 - xxx)
% University of Washington
% EE547 - Final Project
% Fall 2017 - Minseg Robot
% Professor: Professor Linda Bushnell
% TA: Sang UK Sagong

% Step O: Close all, clear all, clc
close all;
clear ;
clc;

% % Part 4.1 Linear Dynamical Model of the Minseg Robot
% 
% fprintf('\n Part 4.1 - Linear Dynamical Model of the Minseg Robot');
% % Step 1: Define the matrices A, B, C, and D to develop a linear continuous
% % time state-space representation of the system.
% 
% fprintf('Step 1: Define the matrices A, B, C, and D to develop a linear \n ');
% fprintf('continuous time state-space representation of the system.\n');
% 
% % Define all the variables
% syms icmw ip L mp mw rw g kb kt R
% dom1 = icmw*(ip+L*L*mp);
% dom2 = (L*L*mp*mw + ip*(mp+mw))*rw*rw;
% 
% % Define the matrices A B C D
% A1 = (g*L*mp*(icmw + (mp+mw)*rw*rw))/(dom1+dom2);
% A2 = -(kb*kt*(icmw+rw*(mw*rw+mp*(L+rw))))/((R*dom1)+dom2);
% A3 = -(kb*kt*(icmw+rw*(mw*rw+mp*(L+rw))))/((R*rw*dom1) + dom2);
% 
% A4 = (g*L*L*mp*mp*rw*rw)/(dom1+dom2);
% A5 = -(kb*kt*rw*(ip+L*mp*(L+rw)))/((R*dom1)+dom2);
% A6 = -(kb*kt*(ip+L*mp*(L+rw)))/((R*dom1)+dom2);
% A = [0 1 0 0;
%     A1 A2 0 A3;
%     0 0 0 1;
%     A4 A5 0 A6];
% 
% B1 = -(kt*(icmw+rw*(mw*rw+mp*(L+rw))))/((R*dom1)+dom2);
% B2 = -(kt*rw*(ip+L*mp*(L+rw)))/((R*dom1)+dom2);
% B = [0;B1;0;B2];
% C = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
% D = [0;0;0;0];

% Step 2: Measure the physical parameters of your Minseg. Create a table to
% list the values of your measurement (SI Unit(. A set reference values of
% (kt, kb, R) are given as (kt,kb,R) = (0.3233 Nm/A, 0.4953 Vs/rad, 5.2628 Ohms)
kt = 0.3233; % Torque constant Nm/A
kb = 0.4953; % Constant Vs/rad
R = 5.2628; % Resistance of motor Ohms
g = 9.806; % Acceleration m/s^2
L = 0.12; % distance between wheel center and reference point over pendulum
%%mwheel = 0.015; % mass of 1 wheel %kg
m2wheels = 0.036; % mass of 2 wheels + axel
mw = m2wheels; % mass of wheel kg (for both wheel)
mbatteries = 0.26; % mass of batteries total %kg
% % mmotor = 117 - mw;
% % mp = mbatteries + mmotor; % mass of pendulum kg
mp=0.382 %mass of pendulem with battery
rw = 0.021; % radius of wheel m
icmw = mw*rw*rw/2; % moment of inertia at center of mass of wheel kg-m^2
ip = mp*L*L; % moment of inertia at reference point of pendulum kg-m^2


% Part 4.1 Linear Dynamical Model of the Minseg Robot

fprintf('\n Part 4.1 - Linear Dynamical Model of the Minseg Robot\n');
% Step 1: Define the matrices A, B, C, and D to develop a linear continuous
% time state-space representation of the system.

fprintf('\nStep 1: Define the matrices A, B, C, and D to develop a linear \n ');
fprintf('continuous time state-space representation of the system.\n');

% Define all the variables
%syms icmw ip L mp mw rw g kb kt R
dom1 = icmw*(ip+L*L*mp);
dom2 = (L*L*mp*mw + ip*(mp+mw))*rw*rw;

% Define the matrices A B C D
A1 = (g*L*mp*(icmw + ((mp+mw)*rw*rw)))/(dom1+dom2);
A2 = -(kb*kt*(icmw+rw*(mw*rw+mp*(L+rw))))/(R*(dom1+dom2));
A3 = -(kb*kt*(icmw+rw*(mw*rw+mp*(L+rw))))/(R*rw*(dom1+dom2));

A4 = (g*L*L*mp*mp*rw*rw)/(dom1+dom2);
A5 = -(kb*kt*rw*(ip+L*mp*(L+rw)))/(R*(dom1+dom2));
A6 = -(kb*kt*(ip+L*mp*(L+rw)))/(R*(dom1+dom2));
A = [0 1 0 0;
    A1 A2 0 A3;
    0 0 0 1;
    A4 A5 0 A6]

B1 = -(kt*(icmw+rw*(mw*rw+mp*(L+rw))))/(R*(dom1+dom2));
B2 = -(kt*rw*(ip+L*mp*(L+rw)))/(R*(dom1+dom2));
B = [0;B1;0;B2]
C = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]
D = zeros(size(C,1),size(B,2));



% Step 3: Find the transfer function matrix of the linearized system
sss = ss(A,B,C,D)
[num,den] = ss2tf(A,B,C,D)

% Step 4: Find the characteristic polynomial and eigenvalues of matrix A
Char = charpoly(A)
eigA = eig(A)

% Step 5: Is the system asymptotically stable? Is it marginally stable?
% Explain why

% if (eigA > 0)
%     fprintf('\n Since some or all the real part of the eigenvalues of matrix A is positive,');
%     fprintf('the system is not asymptotically stable \n');
% else 
%     fprintf('\n Since all real parts of the eigenvalues of matrix A is negative,');
%     fprintf('the system is asymptotically stable \n');
% end

if all(eigA < 0)
    fprintf('\n Since all real parts of the eigenvalues of matrix A is negative,');
    fprintf('the system is asymptotically stable \n');
    
else
    fprintf('\n Since some or all the real part of the eigenvalues of matrix A is positive,'); 
    fprintf('the system is not asymptotically stable \n');

end

% Step 6: Find the poles of the transfer function. Is the system BIBO
% stable? Explain why
poles = eigA;
ttt = [0,0,0,0];
stem(poles, ttt);
title('Plot of the poles of the transfer function');
xlabel('Real parts of poles');
ylabel('Imagine parts');

fprintf('\n Since there are poles/eigenvalues of matrix A in the right hand');
fprintf('\n side of the xy plane, so the system is not BIBO stable\n');

% Part 4.2 Controllability and Observability of the system

fprintf('\n Part 4.2 Controllability and Observability of the system.');

% Step 7: Find the controllability matrix of the lineared system. What is
% the rank of the controllability matrix? Is the linearized system
% completely controllable?
fprintf('\n The controllability matrix of the system is:\n');
%control_matrix = ctrb(sss.a,sss.b)
control_matrix = ctrb(A,B)
fprintf('\n Checking the ranks of the controllability matrix: \n');
rank_control_matrix = rank(control_matrix)
fprintf('\n Checking the length of matrix A: \n');
lengthA = length(A)
if rank(control_matrix) < length(A)
    fprintf('\n Since the ranks of the control matrix is smaller than the size of matrix A,');
    fprintf('the system is not controllable \n');
else 
    fprintf('\n Since the ranks of the control matrix is greater than or equal to the size');
    fprintf('of matrix A, the system is controllable \n');
end

% Step 8: Analyze the observability of the linearized system with the
% output vector as y:=x=[alpha alphadot x xdot]T
fprintf('\n The observability matrix of the system is:');
%obser_matrix = obsv(sss.a,sss.c)
obser_matrix = obsv(A,C)
fprintf('\n Checking the ranks of the observability matrix: \n');
rank_obser_matrix = rank(obser_matrix)
fprintf('\n Checking the length of matrix A: \n');
lengthA = length(A)
if rank(obser_matrix) < length(A)
    fprintf('\n Since the ranks of the observability matrix is smaller than the size of matrix A,');
    fprintf('the system is not observable \n');
else 
    fprintf('\n Since the ranks of the observability matrix is greater than or equal to the');
    fprintf('size of matrix A, the system is observable \n');
end

% Step 9: Transform the linearized system into a controllable canonical
% form and observable canonical form
%Per documentation Matlab by default provides the CCF form from TF 
den_coeff = den(2:end); % deniminator coefficients of H(s)
controlmatrix_inv = [1, den_coeff(1), den_coeff(2), den_coeff(3);
                     0, 1,            den_coeff(1), den_coeff(2);
                     0, 0,            1,            den_coeff(1);
                     0, 0,            0,            1];
Q = control_matrix * controlmatrix_inv;
AA = round(Q\A*Q*1e5)/1e5;
BB = round(Q\B*1e5)/1e5;
CC = round(C*Q*1e5)/1e5;
DD = D;
fprintf('\n Transform the linearized system into a controllable canonical');
fprintf('form and observable canonical form\n');
fprintf('\n The controllable canonical form is:');
controllable_canonical_form = ss(AA,BB,CC,DD) %%% This needs to be w


fprintf('\n The observable canonical form is:');
observable_canonical_form = canon(sss, 'companion')

% Step 10 - Develop a closed-loop state estimator (fulll-domensional
% observer) for the open-loop system (no feedback yet) such that the poles
% of the obserer are stable and that the dynamics of the observer is at
% least 6-8 times faster than the dynamics o the linearized model. Include
% in your report the value of estimator gain, L
fprintf('\n For our project, we chose the dynamics of the observer is at least');
fprintf('6 times faster than the dynamics of the linearized model');
fprintf('\n The poles of the observer is:');
pole_obser = -6*(abs(poles)+1) %take abs of the pole to make sure the new gain is stable and +1 it's not zero
%%%pole_obser = 6*(poles) %mine
fprintf('\n The estimator gain L of the system is:');
gainL = place(transpose(A), transpose(C), pole_obser)'

% Step 11 - Develop a Simulink model of the linearized system (open-loop
% system with full dimensional observer designer above). Add the state
% estimator derived in the previous step to your Simulink model and set the
% initial conditions of the state estimator to x0^ = [0 0 0 0]. Simulate
% the behavior of the system when a unit step u(t) = 1, t >=0 is aplied at
% the input. Plot the estimated state-variables and output variables on the
% same graph.
tspan=0:0.1:10;
% Bsim =  [0 0 0 0 ;B1 0 0 0;0 0 0 0;B2 0 0 0];
xini = [0 0 0 0];
xhatini = [0 0 0 0];
sim('step11',tspan(end));

figure(1)
%subplot(2,1,1)
plot(time,y,'LineWidth',2)
hold on
plot(time,xhat)
hold on
legend('y_1','y_2','y_3','y_4','x_1obv','x_2obv','x_3obv','x_4obv')
title('Estimated state variables and output variables for gain=L & xhatini=[0 0 0 0]');
grid on
% subplot(2,1,2)
% grid on
% plot(time,y,'-.')
% xlabel('time(sec)')
% 
% legend('y_1','y_2','y_3','y_4')

hold off;

% Step 12 - Consider the case when the linearized system is stabilized
% through the use of feedback control. Using the pole placement method,
% develop in Matlab a proportional controller such that the poles of the
% closed-loop system are stable and that the dynamics of the closed-loop
% model is at least 4-6 times faster than the dynamics of the openloop
% model. Include in your report the value of the proportional gain K
fprintf('\n For our project, we chose the dynamics of the closed-loop model is');
fprintf('at least 6 times faster than the dynamics of the open-loop model');
fprintf('\n The poles of the feedback system is:');
poles_feedback = -6*(abs(poles)+1)
fprintf('\n The gain K of the feedback system is:');
gainK = place(A,B,poles_feedback)

% Step 13 - Derive the state-space representation of the closed-loop
% system. Find the characteristic polynomial and the eigenvalues of the
% closed-loop system. Is this closed-loop system asymptotically stable?
ACL = A - B*gainK;
new_sys = ss(ACL,B,C,D);
fprintf('\nThe new characteristic polynomial of the new system is:');
new_charPoly = poly(ACL)
fprintf('\n The new eigenvalues of the new system is:');
new_eigenvalues = eig(ACL)
if all(new_eigenvalues < 0)
    fprintf('\n Since all real parts of the new eigenvalues of the system is negative,');
    fprintf('the system is asymptotically stable \n');
    
else
    fprintf('\n Since some or all the real part of the new eigenvalues of the system is positive,'); 
    fprintf('the system is not asymptotically stable \n');

end


% Step 14 - Develop a Simulink model of the linearized closed-loop system 
% (no estimator here) when the output of the system equals state variables. 
% Add the proportional controller developed above to the Simulink
% model and simulate the response of the closed-loop system when a unit
% step u(t) = 1, t>=0 is applied. Plot all the outputs on the same graph

sim('step14');
figure(2);
plot(t14,y14, 'LineWidth',1.5);
title('Step input responses of the closed-loop system');
xlabel('Time(second)');
ylabel('');
legend ('y_1','y_2','y_3','y_4')

%% Fixing above
% Step 15 - Combine the feedback controller with the state estimator. This
% should be an easy step since both parts are designed separately already.
% Combine the two in Simulink and simulate the system to see how it
% performs. Try using the estimator for states not measured. Plot the error
% function and discuss your result.

% Step 16 - Bonus Step - Demonstrate the MinSeg balancing with one motor.
% This can be done with no batteries (tethered) or with batteries
% (untethered). You can use the given LQR controller or design your own PID
% 











clear all;

clc; 

A = [0      1              0           0;
     0 -0.1195          17.1429        0;
     0      0              0           1;
     0   -0.2915       65.7143         0];
B = [     0;
     1.1953;
          0;
        2.9155];
    
C = [1 0 0 0;
     0 0 1 0];
 
D = [0;
     0];

 
 states = {'x' 'x_dot' 'phi' 'phi_dot'};
inputs = {'u'};
outputs = {'x'; 'phi'};

sys_ss = ss(A,B,C,D,'statename',states,'inputname',inputs,'outputname',outputs)
 
sys_tf = tf(sys_ss);

inputs = {'u'};
outputs = {'x'; 'phi'};

set(sys_tf,'InputName',inputs)
set(sys_tf,'OutputName',outputs)

figure; 
t1=0:0.01:1;
impulse(sys_tf,t1);
title('Open-Loop Impulse Response')

figure;  
t = 0:0.05:10;
u = ones(size(t));
[y,t] = lsim(sys_tf,u,t);
plot(t,y)
title('Open-Loop Step Response')
axis([0 3 0 50])
legend('x','phi')

figure;  
t3 = 0:0.01:1;
u = zeros(size(t3));
[y3,t3] = lsim(sys_tf,u,t3);
plot(t3,y3)
title('Zero Input Response')
axis([0 1 0 50])
legend('x','phi')
clear
close all
% Last updated on: 28. Feb. 2021
% Version of Matlab and Simulink: R2020b


%% 1.Initial conditions
x_initial = [0 0 -27*pi/180 0 0 0] ; 


%% 2.parameters of the helicopter
I_beta = 1.1404 ;
I_gamma = 0.0398 ;
l1 = 0.655;
g= 9.81 ; 
l2 = 0.4905;
e = 0.042;
b = 0.1775;
m1 = 1.5682; 
m_cw = 1.918;


%% 3.Motor force-voltage parameters
p1 = 6.103 ;
p2 = 4.68 ;


%% 4.System
syms u1 u2 alpha beta gamma dalpha dgamma dbeta 

ddalpha = -(u1 + u2)* sin(gamma) * l1 / (I_beta * cos(beta));

ddbeta =((u1+u2)*cos(gamma)*l1- m1*g*(l1*cos(beta)+ e*sin(beta)) + m_cw*g*(l2*cos(beta)- e*sin(beta)))/I_beta ;

ddgamma = (u1-u2)*b/I_gamma ;

dx= [dalpha; ddalpha;  dbeta; ddbeta  ;dgamma ; ddgamma] ;


%% 5.Linearization
x0 = [0 0 -15*pi/180 0 0 0] ;
y_0 =[0 -15*pi/180 0] ;

u0 = zeros(2,1) ;
u0(1) = (m1*g*(l1*cos(x0(3))+e*sin(x0(3)))-m_cw*g*(l2*cos(x0(3))-e*sin(x0(3))))/(2*l1) ;
u0(2) = u0(1) ;

J1 = jacobian(dx,[alpha,dalpha,beta,dbeta,gamma,dgamma]);
J2 = jacobian(dx,[u1,u2]);

A =double(subs(J1,[alpha,beta,gamma,u1,u2],[x0(1),x0(3),x0(5),u0(1),u0(2)]));
B =double(subs(J2,[alpha,beta,gamma],[x0(1),x0(3),x0(5)]));
C = [1 0 0 0 0 0 ;
     0 0 1 0 0 0 ;
     0 0 0 0 1 0 ];
D=zeros(3,2);


%% 6.Controller design
C_new =  zeros(2,6);
C_new(1,1) = 1 ;
C_new(2,3)= 1 ;


% Observer
E=place(A',C', [-23 -24 -25 -26 -27 -28]); 
E=E';

% Weights
R = diag ([1,1]) ;
Q = diag ([1,1,1,1,1,1,1000,1000]) ;
% Q = diag([10,1,10,1,10,1,50,50]);
sys = ss(A,B,C_new,[]);
K = lqi(sys,Q,R);


%% 7.Pre-filter
V = -inv(C_new/(A-B*K(1:2,1:6))*B); 


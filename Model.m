%% Using paper "Control and Estimation of a Quadcopter Dynamical Model"
% transcribing the basic EoM/State-space model

g = 9.87;

%Assuming all MoI are 1 for convenience (this is a big assumption; will be
%changed as more data is known.)
Izbzb = 1;
Iybyb = 1;
Ixbxb = 1;
Ixx = 1;
Iyy = 1;
Izz = 1;
m = 1;

%Assuming all thrusts = 1 for convenience (this is a big assumption, again)
l = 1;

%B equations
b1 = l/Ixx;
b2 = l/Iyy;
b3 = l/Izz;

roll = 0;
roll_rate = 0;
roll_accel_const = 0;

pitch = 0;
pitch_rate = 0;
pitch_accel_const = 0;

yaw = 0;
yaw_rate = 0;
yaw_accel_const = 0;

x = 0;
x_d = 0;

y = 0;
y_d =0;

z = 0;
z_d = 0;

%Engine thrust commands (Engine 1 is top-left, then continue clockwise.)
t1 = 0;
t2 = 0;
t3 = 0;
t4 = 0;

%Dynamics formulae (EoM)
%Theta: Roll, Phi: pitch, Psi: Yaw
roll_accel_formula = (l*((t2+t3)-(t4+t1))-(Izbzb-Iybyb)*pitch_rate*yaw_rate)/Ixbxb;
pitch_accel_formula = (l*((t3+t4)-(t1+t2))-(Ixbxb-Izbzb)*roll_rate*yaw_rate)/Iybyb;
yaw_accel_formula = ((t1+t3)-(t2+t4) - (Iybyb-Ixbxb)*roll_rate*pitch_rate)/Izbzb;

% Matrixes
U = [t1+t2+t3+t4; (t2+t3)-(t4+t1); (t3+t4)-(t1+t2); (t1+t3)-(t2+t4)];
X = [pitch; pitch_rate; roll; roll_rate; yaw; yaw_rate; z; z_d; x; x_d; y; y_d];

%x_d matrix (A)
x_d = zeros(12,12);
x_d(12,1) = -g;
x_d(10,2) = g;
x_d(1,2) = 1;
x_d(3,4) = 1;
x_d(5,6) = 1;
x_d(7,8) = 1;
x_d(9,10) = 1;
x_d(11,12) = 1;

%u_d matrix (B)
u_d = zeros(12,4);
u_d(8,1) = 1/m;
u_d(2,2) = b1;
u_d(4,3) = b2;
u_d(6,4) = b3;

%C matrix
c = zeros(4,12);
c(1,1) = 1;
c(2,3) = 1;
c(3,5) = 1;
c(4,7) = 1;

%D matrix
D = zeros(4,4);

% state-space format
x_dot = x_d*X + u_d*U;
y = c*X;

%sys = ss(x_d,u_d,c,zeros(4,4));

[K,S,e] = lqr(x_d,u_d);
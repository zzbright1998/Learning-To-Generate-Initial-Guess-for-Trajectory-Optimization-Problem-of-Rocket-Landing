clear all
r0 = [1500;0;2000]; %m
V0 = [-75;0;100]; %m/s
m0 = 1905; %kg
Isp = 225; %seconds
g_vec = [-3.7114;0;0];
Thrust_Num = 3100;
T_max = 0.8*Thrust_Num;
T_min = 0.3*Thrust_Num;
no_of_thrusters = 6;
cant_angle = 27*(pi/180);
alpha = 1/(Isp*9.807*cos(cant_angle));

time0 = 0;
final_time = 72; %seconds
time_step = 1;%0.01;
load('UHis.mat');
load('SigmaHis.mat');

%%%% Control Input %%%%%%%%%%%%
time_array = time0:time_step:final_time;
U_input_acc(:,1) = interp1(UHis(:,1),UHis(:,2),time_array,'linear','extrap');
U_input_acc(:,3) = interp1(UHis(:,1),UHis(:,2),time_array,'linear','extrap');
U_input_acc(:,2) = 0*U_input_acc(:,3);
Sigma_input(:,1) = interp1(SigmaHis(:,1),SigmaHis(:,2),time_array,'linear','extrap');
%%%%%%%%%%%%%%%%%%%%
%%%% State Dynamics %%%%%%%%%%%%
StateHis(1,:) = [r0;V0;m0];
k_state = 1;

A_MAT = zeros(7,7);
I_MAT_3 = eye(3);
I_MAT_7 = eye(7);
A_MAT(1:3,4:6) = I_MAT_3;
B_MAT = zeros(7,4);
B_MAT(4:6,1:3) = I_MAT_3;

STM_k = exp(A_MAT*time_step);
STM_INTEG = integral(@(s) expm(A_MAT.*(time_step - s)),0,time_step, 'ArrayValued', true);
STM = 
for time_k = time0:time_step:(final_time-time_step)
    r = StateHis(k_state,1:3);
    V = StateHis(k_state,4:6);
    m = StateHis(k_state,7);
    %control_input_acc(1,:)
    total_acc = interp1(time_array,control_input_acc,time_k,'linear','extrap')';
    total_thrust = 1000*(interp1(time_array,Thrust_input_acc,time_k,'linear','extrap')')';
    
    control_input = total_acc + g_vec;
    
    %r_next = r + time_step*(V);
    %V_next = V + time_step*(total_thrust/m + g_vec');
    %m_next = m + time_step*(-1*alpha*norm(total_thrust));
    %pause%%%%%%%%%%%%%%%%%%%%%%%%%
    r_next = r + time_step*(V);
    V_next = V + time_step*(total_acc)';
    m_next = m + time_step*(-1*alpha*m*norm(total_acc - g_vec));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    theta_deg = acos(dot(total_acc,g_vec));
    %ControlHis(k_state,1:3) =  control_input';
    ControlHis(k_state,1:3) =  m*total_acc;
    k_state = k_state + 1;
    %size(r_dot_vec)
    %size(V_dot_vec)
    StateHis(k_state,1:8) = [r_next, V_next, m_next, theta_deg];
    ThrustHis(k_state,:) = m_next*norm(total_acc - g_vec);
    SigmaHis(k_state,:) = norm(total_acc - g_vec);
end

ControlHis(k_state,1:3) =  control_input';

figure(1);plot(time_array,ControlHis);grid on;hold on;
xlabel('Time (s)','FontSize',12,'FontWeight','normal'); 
ylabel('Acceleration (m/s^2)','FontSize',12,'FontWeight','normal');
legend('Acc in Z', 'acc in Y', 'acc in X');

figure(2);plot(time_array,StateHis(:,1:3));grid on; hold on;
xlabel('Time (s)','FontSize',12,'FontWeight','normal'); 
ylabel('Position (m)','FontSize',12,'FontWeight','normal');
legend('Position in Z', 'Position in Y', 'Position in X');

figure(3);plot(time_array,StateHis(:,7));grid on; hold on;
xlabel('Time (s)','FontSize',12,'FontWeight','normal'); 
ylabel('Mass (kg)','FontSize',12,'FontWeight','normal');

figure(4);plot(time_array,StateHis(:,4:6));grid on; hold on;
xlabel('Time (s)','FontSize',12,'FontWeight','normal'); 
ylabel('Velocity (m/s)','FontSize',12,'FontWeight','normal');
legend('Vel in Z', 'Vel in Y', 'Vel in X');

figure(5);plot(time_array,ThrustHis/(Thrust_Num*no_of_thrusters));grid on; hold on;
xlabel('Time (s)','FontSize',12,'FontWeight','normal'); 
ylabel('Throttle Level','FontSize',12,'FontWeight','normal');
%legend('Vel in Z', 'Vel in Y', 'Vel in X');

figure(6);plot(time_array,SigmaHis);grid on; hold on;
xlabel('Time (s)','FontSize',12,'FontWeight','normal'); 
ylabel('Sigma','FontSize',12,'FontWeight','normal');

% figure(5);plot(time_array,(180/pi)*StateHis(:,8));grid on; hold on;
% %figure(6);plot(time_array,StateHis(:,7).*ControlHis/1000);grid on;hold on;
% figure(6);plot(time_array,ControlHis/1000);grid on;hold on;
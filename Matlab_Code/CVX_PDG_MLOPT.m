clear all
close all
% r0 = [1500;0;2000]; %m
% V0 = [-75;0;100]; %m/s
% Define the ranges for r1, r2, and r3
rng(98);
r1_range = [1250, 1750];
r2_range = [-100, 100];
r3_range = [1750, 2250];

v1_range = [-75,-75];
v2_range = [0,0];
v3_range = [100,100];
MCnum = 6000;

% Draw 600 integer samples for each component
r1_samples = randi(r1_range, MCnum, 1);  % 600 samples between 1250 and 1750
r2_samples = randi(r2_range, MCnum, 1);  % 600 samples between -100 and 100
r3_samples = randi(r3_range, MCnum, 1);  % 600 samples between 1750 and 2250

v1_samples = randi(v1_range, MCnum, 1);  % 600 samples between 1250 and 1750
v2_samples = randi(v2_range, MCnum, 1);  % 600 samples between -100 and 100
v3_samples = randi(v3_range, MCnum, 1);  % 600 samples between 1750 and 2250

% Combine samples into a single 600x3 matrix
r_samples = [r1_samples, r2_samples, r3_samples];
V_samples = [v1_samples, v2_samples, v3_samples];

for MCi=1:1:MCnum
clearvars -except MCi MCnum r_samples V_samples ML_data
MCi
r0 = r_samples(MCi,:)'; %m
V0 = [-75;0;100]; %m/s
%V0(3) = ((-1)^MCi)*V0(3);
final_time = 82;
%%%%%%%%%%%%%%%%%%%%%%%%%%
%V0 = [-75;0;100]
% r0 = [5000;0;0]; %m
% V0 = [0;0;0]; %m/s
%%%%%%%%%%%%%%%%%%%%%%%%%%
global m0 Isp g_vec Thrust_Num T_max T_min no_of_thrusters cant_angle rho_1 rho_2 alpha
m0 = 1905; %kg
Isp = 225; %seconds
g_vec = [-3.7114;0;0];
Thrust_Num = 3100;
T_max = 0.8*Thrust_Num;
T_min = 0.3*Thrust_Num;
no_of_thrusters = 6;
cant_angle = 27*(pi/180);
rho_1 = no_of_thrusters*T_min*cos(cant_angle);
rho_2 = no_of_thrusters*T_max*cos(cant_angle);
alpha = 1/(Isp*9.807*cos(cant_angle));

time0 = 0;
%final_time = 75;%72; %seconds
time_step = 1;%0.01;
time_array = time0:time_step:final_time;
%load('UHis.mat');
%load('SigmaHis.mat');
load('shaped_Eta_time.mat');

%%%% Control Input %%%%%%%%%%%%
% time_array = time0:time_step:final_time;
% time_array_72 = time0:1:72;
% time_array_81 = time0:1:81;
% U_input_acc(:,1) = interp1(time_array_81',UHis(:,1),time_array,'linear','extrap');
% U_input_acc(:,3) = interp1(time_array_81',UHis(:,3),time_array,'linear','extrap');
% U_input_acc(:,2) = U_input_acc(:,2);
% Sigma_input(:,1) = interp1(time_array_81',SigmaHis(:,1),time_array,'linear','extrap');
% Total_Control = [U_input_acc, Sigma_input];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% State Trajectory Guess %%%%%%%%%%%%%%%%%%%%%

N = length(time_array);
guess_param = 1;
[StateHis,Cap_PSI_MAT,XI_MAT,F_XI_Vec,F_PSI_Vec] = StateAndConstraints(r0,V0,m0,time0,time_step,final_time,shaped_Eta_time,guess_param);
% figure(4);plot(time_array(1,:),StateHis(:,1:3));grid on; hold on;
% xlabel('Time (s)','FontSize',12,'FontWeight','normal'); 
% ylabel('Position (m)','FontSize',12,'FontWeight','normal');
% legend('Position in Z', 'Position in Y', 'Position in X');
% 
% pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%% Cost Function %%%%%%%%%%%% (Omega^T)(eta)
%N = length(time_array);
e_Sigma = [0, 0, 0, 1];
e_Sigma_dt = time_step*e_Sigma;

Omega = [];
for i=1:1:N
    Omega = [Omega, e_Sigma_dt];
end
%size(Omega)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%% Norm U < Sigma Constraint %%%%%%%%%%%%

E_U = [eye(3), zeros(3,1)];

for k_norm=1:1:N
    Cap_Gamma_k = zeros(4,4*(N));
    Cap_Gamma_k(:,(k_norm-1)*4+1:4*(k_norm)) = eye(4,4);

    A_SOC = E_U*Cap_Gamma_k;
    b_SOC = zeros(3,1); %% zero
    d_SOC = (e_Sigma*Cap_Gamma_k)';
    gamma_SOC = 0; %% zero
    
    socConstraints_NormU(k_norm) = secondordercone(A_SOC,b_SOC,d_SOC,gamma_SOC);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%% Mu Constraint Upper%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%A_ineq%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% E_Sigma
E_Sigma = zeros(N,4*N);
for k_Sigma = 1:1:N
    E_Sigma(k_Sigma,4*k_Sigma) = 1;
end
%%%
%%% Mu2_MAT
Mu2_MAT = zeros(N,N);
for k_Mu = 1:1:N
    t_k_Mu = time_array(k_Mu);
    z0_k = log(m0-alpha*rho_2*t_k_Mu);
    Mu2_k = rho_2*exp(-z0_k);
    Mu2_MAT(k_Mu,k_Mu) = Mu2_k;
end
%%%
%%% E_Z
E_Z = zeros(N,7*N);
for k_Z = 1:1:N
    E_Z(k_Z,7*k_Z) = 1;
end
%%% E_Z

%%% Cap_PSI
Cap_PSI = Cap_PSI_MAT;
%%% Cap_PSI
Aineq_Upper_Mu = E_Sigma + Mu2_MAT*E_Z*Cap_PSI;
%%%%%%%%%%%%%%%%%%%%%%%B_ineq%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% One_Vec
One_Vec = ones(N,1);
%%% Z_0_Vec
for k_Z = 1:1:N
    t_k_Z = time_array(k_Z);
    z0_k = log(m0-alpha*rho_2*t_k_Z);
    Z_0_Vec(k_Z,1) = z0_k;
end
%%%
%%% XI_MAT
%Already calculated
%%%
Bineq_Upper_Mu = Mu2_MAT*(One_Vec - (E_Z*XI_MAT) + Z_0_Vec);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%% Mu Constraint Lower%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%A_ineq%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for k_Mu_Lo=2:1:N
%     t_k_Mu_Lo = time_array(k_Mu_Lo);
%     z0_k = log(m0-alpha*rho_2*t_k_Mu_Lo);
%     Mu1_k = rho_1*exp(-z0_k);
%     %%%%%%%%%%%
%     e_Sigma = [0, 0, 0, 1];
%     a_Low = 0.5;
%     b_Low = (1 + z0_k);
%     c_Low = 1/(Mu1_k);
%     d_Low = -1-z0_k-0.5*(z0_k^2);
%     F_Vec = [zeros(1,6),1];
%     XI_k = XI_MAT((k_Mu_Lo-1)*7+1:7*(k_Mu_Lo),1);
%     Cap_PSI_k = Cap_PSI_MAT((k_Mu_Lo-1)*7+1:7*(k_Mu_Lo),:);
%     Cap_Gamma_k = zeros(4,4*(N));
%     Cap_Gamma_k(:,(k_Mu_Lo-1)*4+1:4*(k_Mu_Lo)) = eye(4,4);
%     %%%%%%%%%%%
%     Q_Low = a_Low*(F_Vec*Cap_PSI_k)'*(F_Vec*Cap_PSI_k);
%     q_Low = (2*a_Low*(F_Vec*XI_k)*(F_Vec*Cap_PSI_k) - b_Low*F_Vec*Cap_PSI_k - c_Low*e_Sigma*Cap_Gamma_k)';
%     r_Low = -(d_Low + b_Low*F_Vec*XI_k - a_Low*( (F_Vec*XI_k)^2 ) );
%     %%%%%%%%%%%
%     A_SOC = sqrt(a_Low)*(F_Vec*Cap_PSI_k)';
%     %A_SOC = sqrtm(Q_Low);
%     b_SOC = -A_SOC\q_Low;
%     d_SOC = 0; %% zero
%     gamma_SOC = -sqrt(b_SOC'*b_SOC - r_Low);
%     test_var = b_SOC'*b_SOC- r_Low;
%     %socConstraints_Mu_Low(k_Mu_Lo-1) = secondordercone(A_SOC,b_SOC,d_SOC,gamma_SOC);
%     %pause
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for k_Mu_Lo=2:1:N
%     t_k_Mu_Lo = time_array(k_Mu_Lo);
%     z0_k = log(m0-alpha*rho_2*t_k_Mu_Lo);
%     Mu1_k = rho_1*exp(-z0_k);
%     %%%%%%%%%%%
%     e_Sigma = [0, 0, 0, 1];
%     a_Low = 0.5;
%     b_Low = (1 + z0_k);
%     c_Low = 1/(Mu1_k);
%     d_Low = -1-z0_k-0.5*(z0_k^2);
%     F_Vec = [zeros(1,6),1];
%     XI_k = XI_MAT((k_Mu_Lo-1)*7+1:7*(k_Mu_Lo),1);
%     Cap_PSI_k = Cap_PSI_MAT((k_Mu_Lo-1)*7+1:7*(k_Mu_Lo),:);
%     Cap_Gamma_k = zeros(4,4*(N));
%     Cap_Gamma_k(:,(k_Mu_Lo-1)*4+1:4*(k_Mu_Lo)) = eye(4,4);
%     %%%%%%%%%%%
%     A_SOC_Mu = (b_Low\a_Low)*F_Vec*Cap_PSI_k + (c_Low\a_Low)*e_Sigma*Cap_Gamma_k;
%     b_SOC_Mu = -(d_Low\a_Low) - (b_Low\a_Low)*F_Vec*XI_k;
%     d_SOC_Mu = (F_Vec*Cap_PSI_k)'; %% zero
%     gamma_SOC_Mu = -(F_Vec*XI_k);
%     socConstraints_Mu_Low(k_Mu_Lo-1) = secondordercone(A_SOC_Mu,b_SOC_Mu,d_SOC_Mu,gamma_SOC_Mu);
%     %pause
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% E_Sigma
E_Sigma = zeros(N,4*N);
for k_Sigma = 1:1:N
    E_Sigma(k_Sigma,4*k_Sigma) = 1;
end
%%%
%%% Mu1_MAT
Mu1_MAT = zeros(N,N);
for k_Mu = 1:1:N
    t_k_Mu = time_array(k_Mu);
    z0_k = log(m0-alpha*rho_2*t_k_Mu);
    Mu1_k = rho_1*exp(-z0_k);
    Mu1_MAT(k_Mu,k_Mu) = Mu1_k;
end
%%%
%%% E_Z
E_Z = zeros(N,7*N);
for k_Z = 1:1:N
    E_Z(k_Z,7*k_Z) = 1;
end
%%% E_Z

%%% Cap_PSI
Cap_PSI = Cap_PSI_MAT;
%%% Cap_PSI
Aineq_Low_Mu = E_Sigma + Mu1_MAT*E_Z*Cap_PSI;
%%%%%%%%%%%%%%%%%%%%%%%B_ineq%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% One_Vec
One_Vec = ones(N,1);
%%% Z_0_Vec
for k_Z = 1:1:N
    t_k_Z = time_array(k_Z);
    z0_k = log(m0-alpha*rho_2*t_k_Z);
    Z_0_Vec(k_Z,1) = z0_k;
end
%%%
%%% XI_MAT
%Already calculated
%%%
Bineq_Low_Mu = Mu1_MAT*(One_Vec - (E_Z*XI_MAT) + Z_0_Vec);
%% %%%% Z(t) Constraint %%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k_Zt=2:1:N
    t_k_Zt = time_array(k_Zt);
    Cap_Omega_Vec_2_Low(k_Zt-1,1) = log(m0-alpha*rho_2*t_k_Zt);
    Cap_Omega_Vec_1_Upp(k_Zt-1,1) = log(m0-alpha*rho_1*t_k_Zt);
end
%%%%%%%%% Aineq
Aineq_Upp_zt = F_PSI_Vec(2:N,:);
Aineq_Low_zt = -Aineq_Upp_zt;
Aineq_zt = [Aineq_Upp_zt;Aineq_Low_zt];
%%%%%%%%% Bineq
Bineq_Upp_zt = Cap_Omega_Vec_1_Upp - F_XI_Vec(2:N,:);
Bineq_Low_zt = -(Cap_Omega_Vec_2_Low - F_XI_Vec(2:N,:));
Bineq_zt = [Bineq_Upp_zt;Bineq_Low_zt];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%% Glide Constraint %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GS_Theta_max = 86*(pi/180); %% w.r.t Vertical
for k_GS=1:1:N
    t_k_GS = time_array(k_GS);
    XI_k = XI_MAT((k_GS-1)*7+1:7*(k_GS),1);
    Cap_PSI_k = Cap_PSI_MAT((k_GS-1)*7+1:7*(k_GS),:);
    S_MAT = [0, 1, 0, 0, 0, 0, 0;
             0, 0, 1, 0, 0, 0, 0];
    C_Vector = [tan(GS_Theta_max), 0, 0, 0, 0, 0, 0];
    %%%%%%%%%%%
    A_SOC_GS = S_MAT*Cap_PSI_k;
    b_SOC_GS = -S_MAT*XI_k;
    d_SOC_GS = C_Vector*Cap_PSI_k; %% zero
    gamma_SOC_GS = -C_Vector*XI_k;
    socConstraints_GS(k_GS) = secondordercone(A_SOC_GS,b_SOC_GS,d_SOC_GS,gamma_SOC_GS);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%% Terminal Constraints (Equality)%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cap_PSI_N = Cap_PSI_MAT((N-1)*7+1:7*(N),:);
XI_N   = XI_MAT((N-1)*7+1:7*(N),1);
E_Vec = [eye(6,6),zeros(6,1)];
Aeq_Ter = E_Vec*(Cap_PSI_N);
Beq_Ter = zeros(6,1) - E_Vec*XI_N;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%% Solver
A_INEQ_FINAL = [Aineq_Upper_Mu;-Aineq_Low_Mu;Aineq_zt];
B_INEQ_FINAL = [Bineq_Upper_Mu;-Bineq_Low_Mu;Bineq_zt];
%Eta_Optimal = coneprog(Omega',socConstraints_NormU,A_INEQ_FINAL,B_INEQ_FINAL,Aeq_Ter,Beq_Ter);
%Total_SOC = [socConstraints_NormU,socConstraints_GS];
%options = optimoptions('coneprog', 'InitialPoint', 0);
tic
Eta_Optimal = coneprog(Omega',[socConstraints_NormU,socConstraints_GS],A_INEQ_FINAL,B_INEQ_FINAL,Aeq_Ter,Beq_Ter);
toc
%pause
%Eta_Optimal = linprog(Omega',[],[],Aeq_Ter,Beq_Ter);
%pause
%Eta_Optimal = linprog(Omega',A_INEQ_FINAL,B_INEQ_FINAL,Aeq_Ter,Beq_Ter);
clear shaped_Eta shaped_Eta_time StateHis guess_param
shaped_Eta = reshape(Eta_Optimal,[4,length(time_array)]);
shaped_Eta_time (1:4,:) = shaped_Eta;
shaped_Eta_time(5,:) = time_array;
save("shaped_Eta_time.mat","shaped_Eta_time")
%%%%%%%%%%%%%%%%%%
guess_param = 1;
[StateHis,Cap_PSI_MAT,XI_MAT,F_XI_Vec,F_PSI_Vec] = StateAndConstraints(r0,V0,m0,time0,time_step,final_time,shaped_Eta_time,guess_param);

%%%%%%%%%%%%%%%%%%%%



ML_index_data = [shaped_Eta_time',StateHis];
[MLirow,MLicol] = size(ML_index_data);
ML_data(MLirow*(MCi-1)+1:MLirow*MCi,:) = ML_index_data;
save("ML_data.mat","ML_data");

%pause

figure(1);plot(time_array,shaped_Eta(4,:));grid on;hold on;
xlabel('Time (s)','FontSize',12,'FontWeight','normal'); 
ylabel('Sigma','FontSize',12,'FontWeight','normal');
%legend('Acc in Z', 'acc in Y', 'acc in X');

figure(2);plot(time_array(1:end-1),shaped_Eta(4,1:end-1)'.*exp(StateHis(1:end-1,7)));grid on;hold on;
xlabel('Time (s)','FontSize',12,'FontWeight','normal'); 
ylabel('Throttle','FontSize',12,'FontWeight','normal');

figure(3);plot(time_array,shaped_Eta(1:3,:));grid on;hold on;
xlabel('Time (s)','FontSize',12,'FontWeight','normal'); 
ylabel('U','FontSize',12,'FontWeight','normal');
%legend('Acc in Z', 'acc in Y', 'acc in X');


figure(4);plot(StateHis(:,3),StateHis(:,1));grid on;hold on;
xlabel('X (m)','FontSize',12,'FontWeight','normal'); 
ylabel('Z (m)','FontSize',12,'FontWeight','normal');

% figure(1);plot(time_array,ControlHis);grid on;hold on;
% xlabel('Time (s)','FontSize',12,'FontWeight','normal'); 
% ylabel('Acceleration (m/s^2)','FontSize',12,'FontWeight','normal');
% legend('Acc in Z', 'acc in Y', 'acc in X');



figure(5);plot(time_array(1,:),StateHis(:,1:3));grid on; hold on;
xlabel('Time (s)','FontSize',12,'FontWeight','normal'); 
ylabel('Position (m)','FontSize',12,'FontWeight','normal');
legend('Position in Z', 'Position in Y', 'Position in X');

figure(6);plot3(StateHis(:,3),StateHis(:,2),StateHis(:,1));grid on; hold on;
%xlabel('Time (s)','FontSize',12,'FontWeight','normal'); 
%ylabel('Position (m)','FontSize',12,'FontWeight','normal');
%legend('Position in Z', 'Position in Y', 'Position in X');

% figure(3);plot(time_array,StateHis(:,7));grid on; hold on;
% xlabel('Time (s)','FontSize',12,'FontWeight','normal'); 
% ylabel('Mass (kg)','FontSize',12,'FontWeight','normal');
% 
figure(8);plot(time_array,StateHis(:,4:6));grid on; hold on;
xlabel('Time (s)','FontSize',12,'FontWeight','normal'); 
ylabel('Velocity (m/s)','FontSize',12,'FontWeight','normal');
legend('Vel in Z', 'Vel in Y', 'Vel in X');
% 
% figure(5);plot(time_array,ThrustHis/(Thrust_Num*no_of_thrusters));grid on; hold on;
% xlabel('Time (s)','FontSize',12,'FontWeight','normal'); 
% ylabel('Throttle Level','FontSize',12,'FontWeight','normal');
% %legend('Vel in Z', 'Vel in Y', 'Vel in X');
% 
% figure(6);plot(time_array,SigmaHis);grid on; hold on;
% xlabel('Time (s)','FontSize',12,'FontWeight','normal'); 
% ylabel('Sigma','FontSize',12,'FontWeight','normal');

% figure(5);plot(time_array,(180/pi)*StateHis(:,8));grid on; hold on;
% %figure(6);plot(time_array,StateHis(:,7).*ControlHis/1000);grid on;hold on;
% figure(6);plot(time_array,ControlHis/1000);grid on;hold on;

end

writematrix(ML_data, 'ML_data.csv'); 
save("ML_data.mat","ML_data");

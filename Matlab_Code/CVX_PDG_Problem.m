clear all
r0 = [1500;0;2000]; %m
V0 = [-75;0;100]; %m/s
%%%%%%%%%%%%%%%%%%%%%%%%%%
% r0 = [5000;0;0]; %m
% V0 = [0;0;0]; %m/s
%%%%%%%%%%%%%%%%%%%%%%%%%%
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
final_time = 81;%72; %seconds
time_step = 1;%0.01;
%load('UHis.mat');
%load('SigmaHis.mat');
load('shaped_Eta.mat');
UHis = shaped_Eta(1:3,:)';
SigmaHis = shaped_Eta(4,:)';
%%%% Control Input %%%%%%%%%%%%
time_array = time0:time_step:final_time;
time_array_72 = time0:1:72;
time_array_81 = time0:1:81;
% U_input_acc(:,1) = interp1(time_array_72',UHis(:,1),time_array,'linear','extrap');
% U_input_acc(:,3) = interp1(time_array_72',UHis(:,3),time_array,'linear','extrap');
U_input_acc(:,1) = interp1(time_array_81',UHis(:,1),time_array,'linear','extrap');
U_input_acc(:,3) = interp1(time_array_81',UHis(:,3),time_array,'linear','extrap');
U_input_acc(:,2) = U_input_acc(:,2);
%%%%%%%%%%
for i=1:1:length(U_input_acc(:,1))
    thrust_norm(i,1) = norm(U_input_acc(i,:));
end
% Sigma_input(:,1) = interp1(time_array_72',SigmaHis(:,1),time_array,'linear','extrap');
Sigma_input(:,1) = interp1(time_array_81',SigmaHis(:,1),time_array,'linear','extrap');
%pause
Total_Control = [U_input_acc, Sigma_input];

% figure(1);plot(time_array, U_input_acc);
% figure(2);plot(time_array, Sigma_input);
% figure(3);plot(time_array, thrust_norm);
% figure(4);plot(time_array, thrust_norm./Sigma_input);
% pause

%%%%%%%%%%%%%%%%%%%%
%%%% State Dynamics %%%%%%%%%%%%
StateHis(1,:) = [r0;V0;log(m0)];
k_state = 1;
N = length(time_array);
F_Vec = [zeros(1,6),1];
Cap_PSI_MAT((k_state-1)*7+1:7*(k_state),:) = zeros(7,(4*(N)));
XI_MAT((k_state-1)*7+1:7*(k_state),1) = StateHis(1,:);%zeros(7,1);
F_XI_Vec(k_state,1) = F_Vec*XI_MAT((k_state-1)*7+1:7*(k_state),1);
F_PSI_Vec(k_state,:) = F_Vec*Cap_PSI_MAT((k_state-1)*7+1:7*(k_state),:);

A_MAT = zeros(7,7);
I_MAT_3 = eye(3);
I_MAT_7 = eye(7);
A_MAT(1:3,4:6) = I_MAT_3;
B_MAT = zeros(7,4);
B_MAT(4:6,1:3) = I_MAT_3;
B_MAT(7,4) = -alpha;

B_MULT = [];
B_G_MULT = [];
U_MULT = [];
G_MULT = [];
G_Vec = [0;0;0;g_vec;0];

STM_k_A = expm(A_MAT*time_step);
STM_INTEG_B = integral(@(s) expm(A_MAT.*(time_step - s)),0,time_step, 'ArrayValued', true);

%U_input_acc = shaped_Eta(1:3,:)';
%Sigma_input = shaped_Eta(4,:)';
%Total_Control = interp1(time_array_81,shaped_Eta',time_array,'linear','extrap');

for time_k = time0:time_step:(final_time-time_step)
    %r = StateHis(k_state,1:3);
    %V = StateHis(k_state,4:6);
    %m = StateHis(k_state,7);
    
    %control_input_U = interp1(time_array,U_input_acc,time_k,'linear','extrap')';
    %control_input_Sigma = interp1(time_array,Sigma_input,time_k,'linear','extrap')';

    %ControlHis(k_state,1:3) =  control_input_U;
    %ControlHis(k_state,4) =  control_input_Sigma;

    B_MULT = [(STM_k_A)^(k_state-1)*(STM_INTEG_B*B_MAT), B_MULT];
    B_G_MULT = [(STM_k_A)^(k_state-1)*(STM_INTEG_B), B_G_MULT];
    U_MULT = [U_MULT;Total_Control(k_state,:)'];
    %U_MULT = [U_MULT;(Total_Control(k_state,:)' - [g_vec;0])];
    G_MULT = [G_MULT;G_Vec];

    %%%%%%%%%%%%%%%%%%%%%%%%%

    Solution_k =  ((STM_k_A)^(k_state))*(StateHis(1,:)') + ...
                   B_MULT*U_MULT + ...
                   B_G_MULT*G_MULT;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    k_state = k_state + 1;
    StateHis(k_state,1:7) = Solution_k';
    
    %%%%% Constraints Calculation %%%%%%%%%%%%%
    Cap_PSI_MAT((k_state-1)*7+1:7*(k_state),:) = [B_MULT, zeros(7,(4*(N) - 4*(k_state-1) ))];
    XI_MAT((k_state-1)*7+1:7*(k_state),1) = ((STM_k_A)^(k_state))*(StateHis(1,:)') + B_G_MULT*G_MULT;
    %%% z(t) constraint calculation F*XI,F*PSI
    F_XI_Vec(k_state,1) = F_Vec*XI_MAT((k_state-1)*7+1:7*(k_state),1);
    F_PSI_Vec(k_state,:) = F_Vec*Cap_PSI_MAT((k_state-1)*7+1:7*(k_state),:);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%ControlHis(k_state,1:3) =  control_input_U;
%ControlHis(k_state,4) =  control_input_Sigma;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%pause
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
tic
Eta_Optimal = coneprog(Omega',[socConstraints_NormU,socConstraints_GS],A_INEQ_FINAL,B_INEQ_FINAL,Aeq_Ter,Beq_Ter);
toc
%pause
%Eta_Optimal = linprog(Omega',[],[],Aeq_Ter,Beq_Ter);
%pause
%Eta_Optimal = linprog(Omega',A_INEQ_FINAL,B_INEQ_FINAL,Aeq_Ter,Beq_Ter);
shaped_Eta = reshape(Eta_Optimal,[4,length(time_array)]);
%%%%%%%%%%%%%%%%%%

figure(1);plot(time_array,shaped_Eta(4,:));grid on;hold on;
xlabel('Time (s)','FontSize',12,'FontWeight','normal'); 
ylabel('Sigma','FontSize',12,'FontWeight','normal');
%legend('Acc in Z', 'acc in Y', 'acc in X');

figure(11);plot(time_array(1:end-1),shaped_Eta(4,1:end-1)'.*exp(StateHis(1:end-1,7)));grid on;hold on;
xlabel('Time (s)','FontSize',12,'FontWeight','normal'); 
ylabel('Throttle','FontSize',12,'FontWeight','normal');

figure(2);plot(time_array,shaped_Eta(1:3,:));grid on;hold on;
xlabel('Time (s)','FontSize',12,'FontWeight','normal'); 
ylabel('U','FontSize',12,'FontWeight','normal');
%legend('Acc in Z', 'acc in Y', 'acc in X');


figure(3);plot(StateHis(:,3),StateHis(:,1));grid on;hold on;
xlabel('X (m)','FontSize',12,'FontWeight','normal'); 
ylabel('Z (m)','FontSize',12,'FontWeight','normal');

% figure(1);plot(time_array,ControlHis);grid on;hold on;
% xlabel('Time (s)','FontSize',12,'FontWeight','normal'); 
% ylabel('Acceleration (m/s^2)','FontSize',12,'FontWeight','normal');
% legend('Acc in Z', 'acc in Y', 'acc in X');



figure(4);plot(time_array(1,:),StateHis(:,1:3));grid on; hold on;
xlabel('Time (s)','FontSize',12,'FontWeight','normal'); 
ylabel('Position (m)','FontSize',12,'FontWeight','normal');
legend('Position in Z', 'Position in Y', 'Position in X');

figure(5);plot3(StateHis(:,3),StateHis(:,2),StateHis(:,1));grid on; hold on;
%xlabel('Time (s)','FontSize',12,'FontWeight','normal'); 
%ylabel('Position (m)','FontSize',12,'FontWeight','normal');
%legend('Position in Z', 'Position in Y', 'Position in X');

% figure(3);plot(time_array,StateHis(:,7));grid on; hold on;
% xlabel('Time (s)','FontSize',12,'FontWeight','normal'); 
% ylabel('Mass (kg)','FontSize',12,'FontWeight','normal');
% 
% figure(4);plot(time_array,StateHis(:,4:6));grid on; hold on;
% xlabel('Time (s)','FontSize',12,'FontWeight','normal'); 
% ylabel('Velocity (m/s)','FontSize',12,'FontWeight','normal');
% legend('Vel in Z', 'Vel in Y', 'Vel in X');
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
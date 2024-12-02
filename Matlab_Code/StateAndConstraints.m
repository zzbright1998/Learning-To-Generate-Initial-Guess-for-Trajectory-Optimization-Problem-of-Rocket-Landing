%%%%%%%%%%%%%%%%%%%%
function [StateHis,Cap_PSI_MAT,XI_MAT,F_XI_Vec,F_PSI_Vec] = StateAndConstraints(r0,V0,m0,time0,time_step,final_time,shaped_Eta_time,guess_param)

global m0 Isp g_vec Thrust_Num T_max T_min no_of_thrusters cant_angle rho_1 rho_2 alpha

%%%%%%%%%%%%%Control evaluation %%%%%%%%%%%%%

shaped_Eta = shaped_Eta_time(1:4,:);
UHis = shaped_Eta(1:3,:)';
SigmaHis = shaped_Eta(4,:)';
Initial_Guess_time = shaped_Eta_time(5,:);

time_array = time0:time_step:final_time;

U_input_acc(:,1) = interp1(Initial_Guess_time',UHis(:,1),time_array,'linear','extrap');
U_input_acc(:,3) = interp1(Initial_Guess_time',UHis(:,3),time_array,'linear','extrap');
U_input_acc(:,2) = interp1(Initial_Guess_time',UHis(:,2),time_array,'linear','extrap');
Sigma_input(:,1) = interp1(Initial_Guess_time',SigmaHis(:,1),time_array,'linear','extrap');
Total_Control = guess_param*[U_input_acc, Sigma_input];

% figure(2);plot(time_array,shaped_Eta(1:3,:));grid on;hold on;
% xlabel('Time (s)','FontSize',12,'FontWeight','normal'); 
% ylabel('U','FontSize',12,'FontWeight','normal');

%%%% State Dynamics %%%%%%%%%%%%
StateHis(1,:) = [r0;V0;log(m0)];
k_state = 1;
time_array = time0:time_step:final_time;
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

end
%%% IC 1 subsurface travel
%r0 = [1500;0;2000]; %m
%V0 = [-75;0;100]; %m/s

%%% IC 2 No subsurface travel
r0 = [3000;1000;2000]; %m
V0 = [-50;10;100]; %m/s
MC_number = 1;

m0 = 1905; %kg
Isp = 225; %seconds
g_vec = [-3.7114;0;0];
T_max = 0.8*3100;
T_min = 0.3*3100;
no_of_thrusters = 6;
cant_angle = 27*(pi/180);
alpha = 1/(Isp*9.807*cos(cant_angle));

Cd = 1;
Area = 3*2.7;

rf = [0,0,0]; %m
Vf = [0,0,0]; %m/s
af = [0,0,0]; %m/s^2
%%%% State Dynamics %%%%%%%%%%%%
StateHis(1,:) = [r0;V0;m0];
%ControlHis(1,1:3) =  control_input';
%StateHis
%pause
k_state = 1;
a0 = [0.7194, 0, -5.4329];
ControlHis(1,:) = a0;
%ThrustHis(k_state,1:3) = m0*(ControlHis(1,:) + g_vec');
ThrustHis(k_state,1) = norm(m0*(ControlHis(1,:) + g_vec'));

Gamma_0 = 20;
sense_2D = 0;

t_go_0 = Time_to_go_Padhi(r0,V0,a0,rf,Vf,af, sense_2D, Gamma_0);
Subsurface_Travel_Flag = Ground_collision_detection(r0(1),V0(1),a0(1),rf(1),Vf(1),af(1),t_go_0);

if (Subsurface_Travel_Flag == 1)
    t_go_0
    t_go_0 = Ground_collision_avoidance_poly_2(r0(1),V0(1),a0(1),rf(1),Vf(1),af(1),t_go_0);
    %pause
end

time0 = 0;
final_time = t_go_0; %seconds
%pause
time_step = 0.01;

time_array = time0:time_step:final_time;
%DragHis = 0.5*

%pause
for time_k = time0:time_step:(final_time-time_step)
    %%%% Initialization%%%%%%%%%%%
    r = StateHis(k_state,1:3);
    V = StateHis(k_state,4:6);
    m = StateHis(k_state,7);
    total_acc = ControlHis(k_state,:);
    %%% Drag Estimation %%%%%%%%%%%%%%
    T = -31 - 0.000998 * r(1);
    p = 0.699 * exp(-0.00009 * r(1));  % K-Pa
    rho = p / (0.1921 * (T + 273.1)); % kg/m^3
    DragHis(k_state,:) = -1*sign(V).*(0.5*rho*(V.^2)*Cd*Area);
    %%%% State Dynamics Integration%%%%%%%%%%%%%%%%%%%%%%%%%
    r_next = r + time_step*(V);
    V_next = V + time_step*(total_acc + 0*( DragHis(k_state,:)/m ) );
    m_next = m + time_step*(-1*alpha*m*norm(total_acc - g_vec'));
    %%%%% Time-to-go Computation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t_go_0 = final_time - time_k;
    t_go = t_go_0 - time_step;
    %%%%% Acceleration Computation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    total_acc = Guidance_3D_Padhi(r_next, V_next, total_acc, rf, Vf, af, t_go_0, t_go);
    k_state = k_state + 1;
    ControlHis(k_state,1:3) =  total_acc';
    %ThrustHis(k_state,1:3) = m_next*(total_acc + g_vec');
    ThrustHis(k_state,1) = norm(m_next*(total_acc - g_vec'));
    %size(r_dot_vec)
    %size(V_dot_vec)
    StateHis(k_state,1:7) = [r_next, V_next, m_next];
end

DragHis(k_state,:) = DragHis(k_state-1,:);
ControlHis(k_state,1:3) =  total_acc';

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Color_Guess_Matrix = [0.6350 0.0780 0.1840;
                      0.8500 0.3250 0.0980;
                      0.9290 0.6940 0.1250;
                      0.4660 0.6740 0.1880;
                      0.3010 0.7450 0.9330;
                      0 0.4470 0.7410;
                      0.4940 0.1840 0.5560];
                  
Color_Count = rem(MC_number,8);
if (Color_Count == 0)
    Color_Guess = [1 0 1];
else
    Color_Guess = Color_Guess_Matrix(Color_Count,:);
end
Line_St = "-";
Linewidth_Param = 1;
%%

%           % State Constraints %%%%%%%%%%%%%%%
figure(1);plot(time_array,StateHis(:,1), 'Color', Color_Guess,'LineStyle',Line_St,'Linewidth',Linewidth_Param);grid on;hold on;
xlabel('Time (s)','FontSize',12,'FontWeight','normal'); ylabel('Altitude','FontSize',12,'FontWeight','normal');
%figure(1);plot(StateHis(:,2)*(1/1000),7*ones(length(State_Cons_His(:,1)),1), '--r','Linewidth',3);grid on;hold on;
set(gca,'GridAlpha', 0.3);set(gca,'GridColor', 'b');set(gca,'GridLineStyle', '--');
set(gca,'fontname','Times New Roman','FontSize',12);%,'fontangle','italic');
set(figure(1),'renderer','Painters')
%xlim([0.5 6.2]);ylim([0 11]);

figure(2);plot(time_array,StateHis(:,2), 'Color', Color_Guess,'LineStyle',Line_St,'Linewidth',Linewidth_Param);grid on;hold on;
xlabel('Time (s)','FontSize',12,'FontWeight','normal'); ylabel('Crossrange','FontSize',12,'FontWeight','normal');
%figure(2);plot(StateHis(:,2)*(1/1000),8.5*ones(length(State_Cons_His(:,2)),1), '--r','Linewidth',3);grid on;hold on;
set(gca,'GridAlpha', 0.3);set(gca,'GridColor', 'b');set(gca,'GridLineStyle', '--');
%xlim([0.5 6.2]);ylim([0 11]);
set(gca,'fontname','Times New Roman','FontSize',12);%,'fontangle','italic');
set(figure(2),'renderer','Painters')

figure(3);plot(time_array,StateHis(:,3), 'Color', Color_Guess,'LineStyle',Line_St,'Linewidth',Linewidth_Param);grid on;hold on;
xlabel('Time (s)','FontSize',12,'FontWeight','normal'); ylabel('Downrange','FontSize',12,'FontWeight','normal');
%figure(3);plot(StateHis(:,2)*(1/1000),70*ones(length(State_Cons_His(:,2)),1), '--r','Linewidth',3);grid on;hold on;
set(gca,'GridAlpha', 0.3);set(gca,'GridColor', 'b');set(gca,'GridLineStyle', '--');
set(gca,'fontname','Times New Roman','FontSize',12);%,'fontangle','italic');
set(figure(3),'renderer','Painters')

figure(4);plot(time_array,ThrustHis(1:end)/1000, 'Color', Color_Guess,'LineStyle',Line_St,'Linewidth',Linewidth_Param);grid on;hold on;
xlabel('Time (s)','FontSize',12,'FontWeight','normal'); ylabel('Thrust (kN)','FontSize',12,'FontWeight','normal');
%figure(4);plot(TimeLine,-80*ones(length(StateHis(:,7)),1), '--r','Linewidth',3);grid on;hold on;
set(gca,'GridAlpha', 0.3);set(gca,'GridColor', 'b');set(gca,'GridLineStyle', '--');
set(gca,'fontname','Times New Roman','FontSize',12);%,'fontangle','italic');
set(figure(4),'renderer','Painters');
%xlim([0.5 6.2]);ylim([0 11]);
%pause

% figure(1);plot(time_array,ControlHis);grid on;hold on;
% xlabel('Time (s)','FontSize',12,'FontWeight','normal'); ylabel('Acceleration (X,Y,Z)','FontSize',12,'FontWeight','normal');
% figure(2);plot(time_array,StateHis(:,1:3));grid on; hold on;
% xlabel('Time (s)','FontSize',12,'FontWeight','normal'); ylabel('Position (X,Y,Z)','FontSize',12,'FontWeight','normal');
% figure(3);plot(time_array,StateHis(:,7));grid on; hold on;
% xlabel('Time (s)','FontSize',12,'FontWeight','normal'); ylabel('Mass (kg)','FontSize',12,'FontWeight','normal');
% figure(4);plot(time_array,StateHis(:,4:6));grid on; hold on;
% xlabel('Time (s)','FontSize',12,'FontWeight','normal'); ylabel('Velocity (X,Y,Z)','FontSize',12,'FontWeight','normal');
% %figure(5);plot3(StateHis(:,3),StateHis(:,2),StateHis(:,1));grid on; hold on;
% figure(6);plot(time_array(1:end),ThrustHis(1:end)/1000);grid on; hold on;
% xlabel('Time (s)','FontSize',12,'FontWeight','normal'); ylabel('Total Thrust','FontSize',12,'FontWeight','normal');

OutputHis = [ControlHis, StateHis, ThrustHis];
%pause
%figure(7);plot(time_array,DragHis);grid on; hold on;
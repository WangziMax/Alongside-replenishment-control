%% ONRT circle test
% ������ʾEuler�����΢�ַ���
clear;clc;close all
wave_flag=[0 1];
L=3.147;        % ����
ts=0;           % ���濪ʼʱ��
te=120;          % �������ʱ��
dt=0.01;        % ʱ����
t=ts:dt:te;
x=zeros(length(t),6);                       % ״̬�����洢
rudder=0/180*pi;                           % ���
rps=538/60;                                 % ������ת�� rps  ת/��
x(1,:)=[1.1 0 0 0 0 rudder];                % ״̬������ֵ
fw2=[-3, -3, -0.3];                              % 2nd-order wave drift forces (N,N,N*m)
Store_fw2=zeros(length(t),3);
Store_fw2(1,:)=fw2;
psi=zeros(length(t),1);                     % �����
F_wave=zeros(length(t),3);                  % ��������ϵ�µĲ�����
F_wave(1,:)=fw2;
beta_orig   = 45*pi/180;                    % encounter angle in earch coordinates
w0          = 1.1;                          % wave frequency (rad/s)
x_w         = [0 0]';
k1          = 0.03;
k2          = 0.02;
k3          = 0;
for i=2:length(t)
    fw2_dot=[0.1*randn(1) 0.1*randn(1) 0.1*randn(1)];
    fw2= fw2+dt*fw2_dot;
    if wave_flag(2)==0
        x_dot=ONRT(x(i-1,:),rps,[0 0 0]);             % δ�Ӳ�����
    else
        x_dot=ONRT(x(i-1,:),rps,F_wave(i-1,:));
    end
    x(i,1:4)=x(i-1,1:4)+x_dot(1:4)*dt;
    x(i,5)=x(i-1,5)+(x(i-1,3)+x(i,3))*dt/2;
    x(i,6)=rudder;
    psi(i)=psi(i-1)+(x(i-1,4)+x(i,4))*dt/2;             % �����
    F_wave(i,:)=[fw2(1)*cos(psi(i))+fw2(2)*sin(psi(i)) (-fw2(1)*sin(psi(i))+fw2(2)*cos(psi(i)))*0.6 fw2(3)];  % ��������ϵ�µĲ�����
    Store_fw2(i,:)= fw2(:)';
end
% plot
% trajectory
traj_x=zeros(length(t),1);
traj_y=zeros(length(t),1);
for i=2:length(t)
    u_bar=(x(i-1,1)+x(i,1))/2;
    v_bar=(x(i-1,2)+x(i,2))/2;
    if wave_flag(1)==0
        traj_x(i)=traj_x(i-1)+u_bar*dt*cos(psi(i-1))-v_bar*dt*sin(psi(i-1));
        traj_y(i)=traj_y(i-1)+u_bar*dt*sin(psi(i-1))+v_bar*dt*cos(psi(i-1));
    elseif wave_flag(1)==1
        U = sqrt(u_bar^2 + v_bar^2);
        beta = beta_orig- psi(i-1);
        [disw1, x_w] = WaveDisplacement(w0, dt, x_w, beta, U);
        Disw1 = [-k1*cos(beta), -k2*sin(beta), k3]'*disw1;	
        traj_x(i)=traj_x(i-1)+(u_bar*dt+Disw1(1))*cos(psi(i-1))-(v_bar*dt+Disw1(2))*sin(psi(i-1));
        traj_y(i)=traj_y(i-1)+(u_bar*dt+Disw1(1))*sin(psi(i-1))+(v_bar*dt+Disw1(2))*cos(psi(i-1));
    end
end
% ����״̬����ʱ������
figure('units','normalized','position',[0.1 0.3 0.4 0.5])
plot(t,x(:,1),'r',t,x(:,2),'g',t,x(:,4),'b',t,-x(:,5),'c','linewidth',1.2);
grid on
legend('u','v','r','\Phi')
title('state variables')
% �����˶��켣
figure('units','normalized','position',[0.5 0.3 0.4 0.5])
plot(traj_y/L,traj_x/L,'linewidth',1.5);grid on;axis equal
title('trajectory')
xlabel('y')
ylabel('x')


%% ONRT 船模  MMG模型
clear;clc;close all
dbstop if error
ts=0;
te=60;
dt=0.01;
t=ts:dt:te;
x=zeros(length(t),6);
rudder=25;                  % 舵角
x(1,:)=[1.1 0 0 0 0 rudder/180*pi];       % 初始状态变量值
n=538/60;                   % 螺旋桨转速
ship_resistance=4.1656;     % Fr=0.2时的静水直航阻力
iter=1;
for i=2:length(t)
    u=x(i-1,1);
    v=x(i-1,2);
    p=x(i-1,3);
    r=x(i-1,4);
    phi=x(i-1,5);       % 横摇角   rad
    delta=x(i-1,6);     % 舵角     rad
    rho=1000;           % 水密度
    g=9.81;
    U=sqrt(u^2+v^2);    % 合速度
    % 船型参数值
    L=3.147;            B=0.384;            D=0.266;            d=0.112;
    W=72.6;             GM=0.0424;          A_R=0.012*2;        C_b=0.535;
    K_xx=0.444*B;       K_yy=0.246*L;       K_zz=0.246*L;       D_p=0.1066;
    y_G=0.156; % 重心高度（从龙骨算起）
    S_0=1.5;   % 船体湿表面积
    I_x=W*K_xx^2;
    I_z=W*K_zz^2;
    % 水动力系数和舵系数
    epsilong=1.0;                   gamma_R=0.7;                        l_R=-1.0*L;                         t_R=0.3;
    alpha_H=0.25;                   z_HR=0.819*d;                       x_H=-0.45*L;                        m_x=0.01307*1/2*rho*L^2*d;
    X_vv=-0.08577*1/2*rho*L*d;      X_vr=0.052212*1/2*rho*L^2*d;        X_rr=-0.02126*1/2*rho*L^3*d;        m_y=0.10901*1/2*rho*L^2*d;
    Y_v=-0.30015*1/2*rho*U*L*d;     Y_r=-0.08316*1/2*rho*U*L^2*d;       Y_vvv=-1.77272*1/2*rho*L*d/U;       Y_vvr=0.261986*1/2*rho*L^2*d/U;
    Y_vrr=-0.79966*1/2*rho*L^3*d/U; Y_rrr=0.173916*1/2*rho*L^4*d/U;     Y_phi=-0.00051*1/2*rho*U^2*L*d;     J_x=4.1e-5*1/2*rho*L^4*d;
    z_H=0.852*d;                    K_p=-0.2429*1/2*rho*U*L^3*d;        K_phi=0.000626*1/2*rho*U^2*L^2*d;   J_z=0.00789*1/2*rho*L^4*d;
    N_v=-0.09323*1/2*rho*U*L^2*d;   N_r=-0.05494*1/2*rho*U*L^3*d;       N_vvv=-0.53235*1/2*rho*L^2*d/U;     N_vvr=-0.62895*1/2*rho*L^3*d/U;
    N_vrr=-0.13897*1/2*rho*L^4*d/U; N_rrr=-0.00446*1/2*rho*L^5*d/U;     N_phi=-0.00511*1/2*rho*U^2*L^2*d;   t_p=0;
    w_p=0;
    K_temp=z_H*[Y_r,Y_v,Y_rrr,Y_vrr,Y_vvr,Y_vvv];
    K_r=K_temp(1);
    K_v=K_temp(2);
    K_rrr=K_temp(3);
    K_vrr=K_temp(4);
    K_vvr=K_temp(5);
    K_vvv=K_temp(6);
    % 附加
    lambda_R=1.1814;       % 舵的展弦比
    x_R=-0.446*L;       % 舵力作用中心至船体重心的纵向距离
    w_R=0;              % 舵的伴流分数
    % 计算法向舵力
    u_p=(1-w_p)*u;
    elta=2/3;
    kappa=0.6/epsilong;
    f_alpha=6.13*lambda_R/(lambda_R+2.25);
    J=(1-w_p)*u/(n*D_p);
    if J<0.
        disp(num2str(iter));
        warning('u <= 0');
        keyboard;
    end
    K_T_S=0.62702-0.26467*J-0.09665*J^2;
    K_T_P=0.64111-0.27016*J-0.09319*J^2;
    K_T=K_T_S+K_T_P;
    u_R=epsilong*u_p*sqrt(elta*(1+kappa*(sqrt(1+8*K_T/pi/J^2)-1))^2+1-elta);
    v_R=-gamma_R*(v+l_R*r);
    U_R=sqrt(u_R^2+v_R^2);
    alpha_R=delta-gamma_R*v_R/u_R;
    F_N=1/2*rho*A_R*f_alpha*U_R^2*sin(alpha_R);
    % 计算舵分力
    X_R=-(1-t_R)*F_N*sin(delta);
    Y_R=-(1+alpha_H)*F_N*cos(delta)*cos(phi);
    K_R=z_HR*(1+alpha_H)*F_N*cos(delta);
    N_R=-(x_R+alpha_H*x_H)*F_N*cos(delta)*cos(phi);
    % 计算螺旋桨推力
    T=(1-t_p)*rho*n^2*D_p^4*K_T;
    % 计算GZ
    phi_d=abs(phi*180/pi);
    GZ_full=6.061*sin(0.01683*phi_d+1.117)+6.184*sin(0.02267*phi_d+4.152)+0.7097*sin(0.06305*phi_d-0.2555)+...
        0.06261*sin(0.1597*phi_d-7.114)+0.02889*sin(0.2135*phi_d-5.881)+0.004875030028724;
    if abs(phi_d)>180
        disp(num2str(iter));
        error('Phi_d is more than 180 degree');
    end
    GZ=sign(phi)*GZ_full/48.94;
    % 计算状态变量导数
    u_dot=(T-ship_resistance+X_rr*r^2+X_vr*v*r+X_vv*v^2+X_R+(W+m_y)*v*r)/(W+m_x);
    v_dot=(Y_r*r+Y_v*v+Y_phi*phi+Y_rrr*r^3+Y_vrr*r^2*v+Y_vvr*v^2*r+Y_vvv*v^3+Y_R-(W+m_x)*u*r)/(W+m_y);
    p_dot=(m_x*z_H*u*r+K_r*r+K_p*p+K_v*v+K_phi*phi-W*g*GZ+K_rrr*r^3+K_vrr*v*r^2+K_vvr*v^2*r...
        +K_vvv*v^3+K_R)/(I_x+J_x);
    r_dot=(N_r*r+N_v*v+N_phi*phi+N_rrr*r^3+N_vrr*r^2*v+N_vvr*r*v^2+N_vvv*v^3+N_R)/(I_z+J_z);
    % 输出状态变量导数值
    x_dot=[u_dot v_dot p_dot r_dot];
    zhi(i-1,:)=[x_dot X_R Y_R K_R N_R J];
    % 更新
    x(i,1:4)=x(i-1,1:4)+x_dot(1:4)*dt;
    x(i,5)=x(i-1,5)+x(i-1,3)*dt;
    x(i,6)=rudder/180*pi;
    iter=iter+1;
end
%% plot
% trajectory
traj_x=zeros(length(t),1);
traj_y=zeros(length(t),1);
psi=zeros(length(t),1);
for i=2:length(t)
    psi(i)=psi(i-1)+x(i-1,4)*dt;             % 首向角
    u_bar=x(i-1,1);
    v_bar=x(i-1,2);
    traj_x(i)=traj_x(i-1)+u_bar*dt*cos(psi(i-1))-v_bar*dt*sin(psi(i-1));
    traj_y(i)=traj_y(i-1)+u_bar*dt*sin(psi(i-1))+v_bar*dt*cos(psi(i-1));
end
% 绘制状态变量时历曲线
figure('units','normalized','position',[0.1 0.3 0.4 0.5])
plot(t,x(:,1),'r',t,x(:,2),'g',t,x(:,4),'b',t,x(:,5),'c:','linewidth',1.2);
grid on
legend('u','v','r','\Phi')
title('state variables')
% 绘制运动轨迹
figure('units','normalized','position',[0.5 0.3 0.4 0.5])
plot(traj_x/L,traj_y/L,'linewidth',1.5);grid on;axis equal;
xlabel('x/L');ylabel('y/L');
title('trajectory')
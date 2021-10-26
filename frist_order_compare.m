%% 两种方法的对比仿真
clc;clear all;close all;
n=2000000;  
step=0.00001;
%% 初始条件 
x2 = 1; 
x3 = 1; 
beta20=2.1;
beta21=0.1;
theta2_hat=0;
theta3_hat=0;
rho2_hat=0.1;
rho3_hat=0.1;
%% 设计参数   
k2=1;  
k3=1;    
r2=0.1;
r21=0.1;
r3=0.1;
r31=0.1;
delta2=0.1;
delta3=0.1;




%% 循环-收敛到残差集
tt=0:step:step*n;   
for i=1:1:n
    t=i*step;
    %% 时变参数 
    bt2=   2+1*sign(sin(2*t));
    thetat2=2+1*sign(sin(1*t));
    
    beta2=(beta20-beta21)*exp(-0.4*t)+beta21;
    dbeta2=-0.4*(beta20-beta21)*exp(-0.4*t);
    %%
    Pi2=pi/(2*beta2)*(1+(tan(pi*x2/2/beta2))^2);
    z21=tan(pi*x2/2/beta2);    
%% 控制
W2=2*beta2/pi/z21*atan(z21);
kappa2= (k2/Pi2+0.5*delta2+1+0.5*W2^2*((theta2_hat-dbeta2/beta2)^2+delta2));
bar_u2=-kappa2*z21;
u2=rho2_hat*bar_u2; 
%% 自适应 1
dtheta2_hat=r2*Pi2*z21*x2;
theta2_hat =  theta2_hat + dtheta2_hat * step; 
%% 自适应 2
drho_hat2=-r21*Pi2*z21*bar_u2;
rho2_hat = rho2_hat + drho_hat2 * step;  
%% 系统模型
    dx2 = bt2*u2+thetat2*x2;
    x2 = x2 + dx2 * step;      
%% 参数采集
    x2s(i)=x2;
    u2s(i) = u2;
    bs(i)=beta2;
    theta2s(i)=theta2_hat;
    rho2s(i)=rho2_hat;    
    bs(i)=beta2;
    bs1(i)=-beta2;
    bbs(i)=bt2;
    the1s(i)=thetat2;
end

%% 循环-TAC主编方法
tt=0:step:step*n;   
for i=1:1:n
    t=i*step; 
    %% 时变参数
     bt3=2+1*sign(sin(2*t));
    thetat2=2+1*sign(sin(1*t));
%% 控制 
kappa3= (k3 +  delta3 +0.5 *(theta3_hat^2+1)) *x3;
bar_u3=-kappa3*x3 ;
u3=rho3_hat*bar_u3;
%% 自适应 1
dtheta3_hat=r3* x3^2;
theta3_hat =  theta3_hat + dtheta3_hat * step; 
%% 自适应 2
drho3_hat=-r31 *x3*bar_u3;
rho3_hat = rho3_hat + drho3_hat * step;  
%% 系统模型
    dx3 = bt3*u3+thetat2*x3;
    x3 = x3 + dx3 * step;  
%% 参数采集
    x3s(i)=x3;
    u3s(i) = u3;
    ts(i)=t; 
    theta3s(i)=theta3_hat;
    rho3s(i)=rho3_hat; 
    bbbs(i)=bt3;
    thes(i)=thetat2;
end

 
figure(1)   
subplot(2,2,1)
plot(  ts,x3s,ts,x2s,ts,bs,ts,bs1,'linewidth',1); 
legend({ 'Controller $1^{[1]}$','The Proposed Controller'  },'interpreter','latex');
ylabel('$x(t)$','interpreter','latex','linewidth',18);  
xlabel('Time $ (s)$','interpreter','latex','linewidth',18);  
grid on;   
subplot(2,2,2)
plot(  ts,u3s,ts,u2s, 'linewidth',1); 
legend({'Controller $1^{[1]}$','The Proposed Controller' },'interpreter','latex');
ylabel('$u(t)$','interpreter','latex','linewidth',18); 
xlabel('Time $ (s)$','interpreter','latex','linewidth',18);  
grid on;  
subplot(2,2,3)
plot( ts,theta3s,ts,rho3s,ts,theta2s,ts,rho2s, 'linewidth',1); 
legend({'$\hat{\theta}$ in Controller $1^{[1]}$','$\hat{\rho}$ in Controller $1^{[1]}$','$\hat{\theta}$ in  the proposed Controller', '$\hat{\rho}$ in   the proposed Controller' },'interpreter','latex');
ylabel('Adaptive Parameters','interpreter','latex','linewidth',18); 
xlabel('Time $ (s)$','interpreter','latex','linewidth',18);  
grid on; 
subplot(2,2,4)
plot( ts,thes,ts,bbbs, 'linewidth',1); 
legend({'$\theta(t)$ ','$b(t)$ ' },'interpreter','latex');
ylabel('System Parameters','interpreter','latex','linewidth',18); 
xlabel('Time $ (s)$','interpreter','latex','linewidth',18);  
grid on; 

 
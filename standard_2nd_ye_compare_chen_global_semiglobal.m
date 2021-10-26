%% 二阶系统仿真-时变自适应-chen方法-YE方法-全局半全局-xxxxxxxxxxxxxxxxxxx
clc;clear all;close all;
n=10000;  
step=0.001;
%%   共用参数--------------------------------------
Arhohat = 0.25 ; %自适应参数初始化
Brhohat = 0.25 ;
Crhohat = 0.25 ;
Athetahat = 0 ; %自适应参数初始化  
Bthetahat =0 ;
Cthetahat = 0 ; %自适应参数初始化   
k1 = 0.1;   %反馈增益
k2 = 0.1;   %反馈增益 
gamma = 0.1  ; %学习系数 
delta = 1  ;  %无需知道delta具体值 
L = 0.4 ;   %衰减系数
beta1 = 0.1;  %最终收敛半径
Ax1 =  1; 
Ax2 = -1; 
Bx1 =  1; 
Bx2 = -1;
Cx1 =  1; 
Cx2 = -1; 
%% 初始条件/设计参数--------------chen-------------------chen----------------------chen-----------
tt=0:step:step*n;   
for i=1:1:n
    t=i*step; 
    Abt = 2  + 2*sign(Ax1*Ax2);
    Athetat = 2  + sin(Ax1*Ax2) +0.8*sin(1*t)   +0.2*sin(Ax1*t) + sign(sin(t));  
    Az1 = Ax1;      
    AW1= 1;
    Akappa1 = (k1+Athetahat*AW1+AW1+1);  
    Aalpha1 = -Akappa1 * Az1;    
    Adalpha1_dx1 = - k1 - Athetahat - 2;    
    Adalpha1_dthetahat = -Ax1; 
    AW2 = [-Adalpha1_dx1;    0];   
    Apsibar = [1+ Adalpha1_dx1*Akappa1 - Adalpha1_dthetahat*gamma *Ax1 - (Adalpha1_dx1*Athetahat);  
              -Adalpha1_dx1+Adalpha1_dthetahat*gamma*Adalpha1_dx1*Ax1];  
    Az2 = Ax2 - Aalpha1;
    Akappa2 = k2+0.5*( delta * (norm(AW2,2))+ delta  + 1 + (norm(Apsibar,2))  );  
    Aubar = -Akappa2 * Az2;
    Au = Arhohat * Aubar; 
    Adthetahat = gamma   * Ax1 * Az1  - gamma * Adalpha1_dx1 * Ax1 *Az2;
    Athetahat = Athetahat + Adthetahat * step; 
    Adrhohat = -gamma  * Az2 * Aubar;
    Arhohat = Arhohat + Adrhohat * step;  
    Adx1 = Ax2+Athetat*Ax1;
    Adx2 = Abt*Au;
    Ax1 = Ax1 + Adx1 * step;    
    Ax2 = Ax2 + Adx2 * step;
    % 参数采集了
    ts(i)=t; 
    Ax1s(i)=Ax1;
    Ax2s(i)=Ax2;
    Aus(i) = Au;  
    Athetas(i)=Athetahat;
    Arhos(i)=Arhohat;
    Athetats(i)=Athetat;
    Abts(i)=Abt;
end


%% ----------------------------------------------YE-semiglobal----------------------------
Bbeta0 = 4.1;   
tt=0:step:step*n;   
for i=1:1:n
    t=i*step;
    Bbeta = (Bbeta0 - beta1) * exp(-L*t) + beta1; 
    Bdbeta = -L * Bbeta0 * exp(-L*t); 
    Bbt = 2 + 0.1 *cos(Bx1)+ sign(Bx1*Bx2);
    Bthetat = 2   + sin(Bx1*Bx2) +0.8*sin(1*t) +0.2*sin(Bx1*t)   + sign(sin(t));   
    Bz1 = tan(pi*Bx1/2/Bbeta);    
    BPi = pi/(2*Bbeta)*(1+Bz1^2);  
    BW1= Bx1/ ( ((pi*Bx1)/(2*Bbeta))  );
    Bkappa1 = k1 / BPi + delta/2*BW1^2 + delta/2 + delta/2/BPi;
    Balpha1 = -Bkappa1 * Bz1 + Bx1*Bdbeta /Bbeta -  Bx1 * Bthetahat;    
    Bdalpha1_dx1 =tan((Bx1*pi)/(2*(beta1 + exp(-L*t)*(Bbeta0 - beta1))))*((delta*tan((Bx1*pi)/(2*(beta1 + exp(-L*t)*(Bbeta0 - beta1))))*(2*beta1 + 2*exp(-L*t)*(Bbeta0 - beta1)))/(2*(beta1 + exp(-L*t)*(Bbeta0 - beta1))*(tan((Bx1*pi)/(2*(beta1 + exp(-L*t)*(Bbeta0 - beta1))))^2 + 1)) + (k1*tan((Bx1*pi)/(2*(beta1 + exp(-L*t)*(Bbeta0 - beta1))))*(2*beta1 + 2*exp(-L*t)*(Bbeta0 - beta1)))/((beta1 + exp(-L*t)*(Bbeta0 - beta1))*(tan((Bx1*pi)/(2*(beta1 + exp(-L*t)*(Bbeta0 - beta1))))^2 + 1))) - Bthetahat - (pi*(tan((Bx1*pi)/(2*(beta1 + exp(-L*t)*(Bbeta0 - beta1))))^2 + 1)*(delta/2 + (delta*(2*beta1 + 2*exp(-L*t)*(Bbeta0 - beta1))^2)/(2*pi^2) + (delta*(2*beta1 + 2*exp(-L*t)*(Bbeta0 - beta1)))/(2*pi*(tan((Bx1*pi)/(2*(beta1 + exp(-L*t)*(Bbeta0 - beta1))))^2 + 1)) + (k1*(2*beta1 + 2*exp(-L*t)*(Bbeta0 - beta1)))/(pi*(tan((Bx1*pi)/(2*(beta1 + exp(-L*t)*(Bbeta0 - beta1))))^2 + 1))))/(2*(beta1 + exp(-L*t)*(Bbeta0 - beta1))) - (Bbeta0*L*exp(-L*t))/(beta1 + exp(-L*t)*(Bbeta0 - beta1));    
    Bdalpha1_dbeta = (Bx1*pi*(tan((pi*Bx1)/(2*Bbeta))^2 + 1)*(delta/2 + (2*Bbeta^2*delta)/pi^2 + (Bbeta*delta)/(pi*(tan((pi*Bx1)/(2*Bbeta))^2 + 1)) + (2*Bbeta*k1)/(pi*(tan((pi*Bx1)/(2*Bbeta))^2 + 1))))/(2*Bbeta^2) - (Bdbeta*Bx1)/Bbeta^2 - tan((pi*Bx1)/(2*Bbeta))*(delta/(pi*(tan((pi*Bx1)/(2*Bbeta))^2 + 1)) + (2*k1)/(pi*(tan((pi*Bx1)/(2*Bbeta))^2 + 1)) + (4*Bbeta*delta)/pi^2 + (Bx1*delta*tan((pi*Bx1)/(2*Bbeta)))/(Bbeta*(tan((pi*Bx1)/(2*Bbeta))^2 + 1)) + (2*Bx1*k1*tan((pi*Bx1)/(2*Bbeta)))/(Bbeta*(tan((pi*Bx1)/(2*Bbeta))^2 + 1)));
    Bdalpha1_dthetahat = -Bx1; 
    BW2 = [Bdalpha1_dx1*BW1;
                         0];
    Bpsibar = [Bdalpha1_dbeta * Bdbeta /Bz1 +  Bdalpha1_dx1*(-Bkappa1-BW1*Bthetahat) + Bdalpha1_dthetahat*gamma*BPi*Bx1+BPi+Bdalpha1_dx1*BW1*Bthetahat; 
              Bdalpha1_dx1-Bdalpha1_dthetahat*gamma*Bdalpha1_dx1];       
    Bz2 = Bx2 - Balpha1;
    Bkappa2 = k2+0.5*( delta * (norm(BW2,2))  + delta + 1 + (norm(Bpsibar,2))  );  
    Bubar = -Bkappa2 * Bz2;
    Bu = Brhohat * Bubar; 
    Bdthetahat = gamma * BPi * Bx1 * Bz1  - gamma * Bdalpha1_dx1 * Bx1 *Bz2;
    Bthetahat = Bthetahat + Bdthetahat * step; 
    Bdrhohat = -gamma  * Bz2 * Bubar;
    Brhohat = Brhohat + Bdrhohat * step;  
    Bdx1 = Bx2+Bthetat*Bx1;
    Bdx2 = Bbt*Bu;
    Bx1 = Bx1 + Bdx1 * step;    
    Bx2 = Bx2 + Bdx2 * step;   
    Bx1s(i)=Bx1;
    Bx2s(i)=Bx2;
    Bus(i) = Bu;
    Bbs(i)=Bbeta;
    Bbs1(i)=-Bbeta;   
    Bthetas(i)=Bthetahat;
    Brhos(i)=Brhohat;
    Bthetats(i)=Bthetat;
    Bbts(i)=Bbt;
end


%% -------------------------------------------YE-------global-------------------------------

Cbeta0 = 1;    
tt=0:step:step*n;   
for i=1:1:n
    t=i*step;
    Cbeta = (Cbeta0 - beta1) * exp(-L*t) + beta1; 
    Cdbeta = -L * Cbeta0 * exp(-L*t); 
    Cboundary=atanh(Cbeta); 
    Cbt = 2 + 0.1 * cos(Cx1)+ sign(Cx1*Cx2);
    Cthetat = 2  + sin(Cx1*Cx2) +0.8*sin(1*t) +0.2*sin(Cx1*t) +sign(sin(t));   
    Cz1 = ( Cbeta * tanh(Cx1))/(Cbeta^2-(tanh(Cx1))^2 );    
    CPi = ( (Cbeta^2-(tanh(Cx1))^2)*Cbeta*(sech(Cx1))^2+2*Cbeta*(tanh(Cx1))^2*(sech(Cx1))^2 )/( (Cbeta^2-(tanh(Cx1))^2)^2 ); 
    CPsi = ( Cdbeta*tanh(Cx1)-2*Cbeta^2*Cdbeta*tanh(Cx1) )/( (Cbeta^2-(tanh(Cx1))^2)^2 );
    CW1 = Cx1/Cz1;   
    Ckappa1 = k1 / CPi + Cthetahat*CW1+CPsi*Cz1/CPi+CW1+1/CPi;
    Calpha1 = -Ckappa1 * Cz1 ;   
    Cdalpha1_dx1 = (2*tanh(Cx1)^2*(tanh(Cx1)^2 - 1)*((9*exp(-(2*t)/5))/10 + 1/10)*((tanh(Cx1)^2 - ((9*exp(-(2*t)/5))/10 + 1/10)^2)^2/((tanh(Cx1)^2*((9*exp(-(2*t)/5))/5 + 1/5))/cosh(Cx1)^2 - ((tanh(Cx1)^2 - ((9*exp(-(2*t)/5))/10 + 1/10)^2)*((9*exp(-(2*t)/5))/10 + 1/10))/cosh(Cx1)^2) + (k1*(tanh(Cx1)^2 - ((9*exp(-(2*t)/5))/10 + 1/10)^2)^2)/((tanh(Cx1)^2*((9*exp(-(2*t)/5))/5 + 1/5))/cosh(Cx1)^2 - ((tanh(Cx1)^2 - ((9*exp(-(2*t)/5))/10 + 1/10)^2)*((9*exp(-(2*t)/5))/10 + 1/10))/cosh(Cx1)^2) - (Cx1*(tanh(Cx1)^2 - ((9*exp(-(2*t)/5))/10 + 1/10)^2))/(tanh(Cx1)*((9*exp(-(2*t)/5))/10 + 1/10)) - (Cthetahat*Cx1*(tanh(Cx1)^2 - ((9*exp(-(2*t)/5))/10 + 1/10)^2))/(tanh(Cx1)*((9*exp(-(2*t)/5))/10 + 1/10)) + (tanh(Cx1)*((2*exp(-(2*t)/5)*tanh(Cx1))/5 - (4*exp(-(2*t)/5)*tanh(Cx1)*((9*exp(-(2*t)/5))/10 + 1/10)^2)/5)*((9*exp(-(2*t)/5))/10 + 1/10))/(((tanh(Cx1)^2*((9*exp(-(2*t)/5))/5 + 1/5))/cosh(Cx1)^2 - ((tanh(Cx1)^2 - ((9*exp(-(2*t)/5))/10 + 1/10)^2)*((9*exp(-(2*t)/5))/10 + 1/10))/cosh(Cx1)^2)*(tanh(Cx1)^2 - ((9*exp(-(2*t)/5))/10 + 1/10)^2))))/(tanh(Cx1)^2 - ((9*exp(-(2*t)/5))/10 + 1/10)^2)^2 - ((tanh(Cx1)^2 - 1)*((9*exp(-(2*t)/5))/10 + 1/10)*((tanh(Cx1)^2 - ((9*exp(-(2*t)/5))/10 + 1/10)^2)^2/((tanh(Cx1)^2*((9*exp(-(2*t)/5))/5 + 1/5))/cosh(Cx1)^2 - ((tanh(Cx1)^2 - ((9*exp(-(2*t)/5))/10 + 1/10)^2)*((9*exp(-(2*t)/5))/10 + 1/10))/cosh(Cx1)^2) + (k1*(tanh(Cx1)^2 - ((9*exp(-(2*t)/5))/10 + 1/10)^2)^2)/((tanh(Cx1)^2*((9*exp(-(2*t)/5))/5 + 1/5))/cosh(Cx1)^2 - ((tanh(Cx1)^2 - ((9*exp(-(2*t)/5))/10 + 1/10)^2)*((9*exp(-(2*t)/5))/10 + 1/10))/cosh(Cx1)^2) - (Cx1*(tanh(Cx1)^2 - ((9*exp(-(2*t)/5))/10 + 1/10)^2))/(tanh(Cx1)*((9*exp(-(2*t)/5))/10 + 1/10)) - (Cthetahat*Cx1*(tanh(Cx1)^2 - ((9*exp(-(2*t)/5))/10 + 1/10)^2))/(tanh(Cx1)*((9*exp(-(2*t)/5))/10 + 1/10)) + (tanh(Cx1)*((2*exp(-(2*t)/5)*tanh(Cx1))/5 - (4*exp(-(2*t)/5)*tanh(Cx1)*((9*exp(-(2*t)/5))/10 + 1/10)^2)/5)*((9*exp(-(2*t)/5))/10 + 1/10))/(((tanh(Cx1)^2*((9*exp(-(2*t)/5))/5 + 1/5))/cosh(Cx1)^2 - ((tanh(Cx1)^2 - ((9*exp(-(2*t)/5))/10 + 1/10)^2)*((9*exp(-(2*t)/5))/10 + 1/10))/cosh(Cx1)^2)*(tanh(Cx1)^2 - ((9*exp(-(2*t)/5))/10 + 1/10)^2))))/(tanh(Cx1)^2 - ((9*exp(-(2*t)/5))/10 + 1/10)^2) - (tanh(Cx1)*((9*exp(-(2*t)/5))/10 + 1/10)*((tanh(Cx1)^2 - ((9*exp(-(2*t)/5))/10 + 1/10)^2)/(tanh(Cx1)*((9*exp(-(2*t)/5))/10 + 1/10)) - ((tanh(Cx1)^2 - ((9*exp(-(2*t)/5))/10 + 1/10)^2)^2*((2*tanh(Cx1)*(tanh(Cx1)^2 - 1)*((9*exp(-(2*t)/5))/5 + 1/5))/cosh(Cx1)^2 - (2*tanh(Cx1)*(tanh(Cx1)^2 - 1)*((9*exp(-(2*t)/5))/10 + 1/10))/cosh(Cx1)^2 + (2*sinh(Cx1)*tanh(Cx1)^2*((9*exp(-(2*t)/5))/5 + 1/5))/cosh(Cx1)^3 - (2*sinh(Cx1)*(tanh(Cx1)^2 - ((9*exp(-(2*t)/5))/10 + 1/10)^2)*((9*exp(-(2*t)/5))/10 + 1/10))/cosh(Cx1)^3))/((tanh(Cx1)^2*((9*exp(-(2*t)/5))/5 + 1/5))/cosh(Cx1)^2 - ((tanh(Cx1)^2 - ((9*exp(-(2*t)/5))/10 + 1/10)^2)*((9*exp(-(2*t)/5))/10 + 1/10))/cosh(Cx1)^2)^2 - (2*Cx1*(tanh(Cx1)^2 - 1))/((9*exp(-(2*t)/5))/10 + 1/10) - (2*Cthetahat*Cx1*(tanh(Cx1)^2 - 1))/((9*exp(-(2*t)/5))/10 + 1/10) - (k1*(tanh(Cx1)^2 - ((9*exp(-(2*t)/5))/10 + 1/10)^2)^2*((2*tanh(Cx1)*(tanh(Cx1)^2 - 1)*((9*exp(-(2*t)/5))/5 + 1/5))/cosh(Cx1)^2 - (2*tanh(Cx1)*(tanh(Cx1)^2 - 1)*((9*exp(-(2*t)/5))/10 + 1/10))/cosh(Cx1)^2 + (2*sinh(Cx1)*tanh(Cx1)^2*((9*exp(-(2*t)/5))/5 + 1/5))/cosh(Cx1)^3 - (2*sinh(Cx1)*(tanh(Cx1)^2 - ((9*exp(-(2*t)/5))/10 + 1/10)^2)*((9*exp(-(2*t)/5))/10 + 1/10))/cosh(Cx1)^3))/((tanh(Cx1)^2*((9*exp(-(2*t)/5))/5 + 1/5))/cosh(Cx1)^2 - ((tanh(Cx1)^2 - ((9*exp(-(2*t)/5))/10 + 1/10)^2)*((9*exp(-(2*t)/5))/10 + 1/10))/cosh(Cx1)^2)^2 + (Cthetahat*(tanh(Cx1)^2 - ((9*exp(-(2*t)/5))/10 + 1/10)^2))/(tanh(Cx1)*((9*exp(-(2*t)/5))/10 + 1/10)) + (4*tanh(Cx1)*(tanh(Cx1)^2 - ((9*exp(-(2*t)/5))/10 + 1/10)^2)*(tanh(Cx1)^2 - 1))/((tanh(Cx1)^2*((9*exp(-(2*t)/5))/5 + 1/5))/cosh(Cx1)^2 - ((tanh(Cx1)^2 - ((9*exp(-(2*t)/5))/10 + 1/10)^2)*((9*exp(-(2*t)/5))/10 + 1/10))/cosh(Cx1)^2) + ((tanh(Cx1)^2 - 1)*((2*exp(-(2*t)/5)*tanh(Cx1))/5 - (4*exp(-(2*t)/5)*tanh(Cx1)*((9*exp(-(2*t)/5))/10 + 1/10)^2)/5)*((9*exp(-(2*t)/5))/10 + 1/10))/(((tanh(Cx1)^2*((9*exp(-(2*t)/5))/5 + 1/5))/cosh(Cx1)^2 - ((tanh(Cx1)^2 - ((9*exp(-(2*t)/5))/10 + 1/10)^2)*((9*exp(-(2*t)/5))/10 + 1/10))/cosh(Cx1)^2)*(tanh(Cx1)^2 - ((9*exp(-(2*t)/5))/10 + 1/10)^2)) + (4*k1*tanh(Cx1)*(tanh(Cx1)^2 - ((9*exp(-(2*t)/5))/10 + 1/10)^2)*(tanh(Cx1)^2 - 1))/((tanh(Cx1)^2*((9*exp(-(2*t)/5))/5 + 1/5))/cosh(Cx1)^2 - ((tanh(Cx1)^2 - ((9*exp(-(2*t)/5))/10 + 1/10)^2)*((9*exp(-(2*t)/5))/10 + 1/10))/cosh(Cx1)^2) + (Cx1*(tanh(Cx1)^2 - ((9*exp(-(2*t)/5))/10 + 1/10)^2)*(tanh(Cx1)^2 - 1))/(tanh(Cx1)^2*((9*exp(-(2*t)/5))/10 + 1/10)) + (tanh(Cx1)*((2*exp(-(2*t)/5)*(tanh(Cx1)^2 - 1))/5 - (4*exp(-(2*t)/5)*(tanh(Cx1)^2 - 1)*((9*exp(-(2*t)/5))/10 + 1/10)^2)/5)*((9*exp(-(2*t)/5))/10 + 1/10))/(((tanh(Cx1)^2*((9*exp(-(2*t)/5))/5 + 1/5))/cosh(Cx1)^2 - ((tanh(Cx1)^2 - ((9*exp(-(2*t)/5))/10 + 1/10)^2)*((9*exp(-(2*t)/5))/10 + 1/10))/cosh(Cx1)^2)*(tanh(Cx1)^2 - ((9*exp(-(2*t)/5))/10 + 1/10)^2)) - (2*tanh(Cx1)^2*(tanh(Cx1)^2 - 1)*((2*exp(-(2*t)/5)*tanh(Cx1))/5 - (4*exp(-(2*t)/5)*tanh(Cx1)*((9*exp(-(2*t)/5))/10 + 1/10)^2)/5)*((9*exp(-(2*t)/5))/10 + 1/10))/(((tanh(Cx1)^2*((9*exp(-(2*t)/5))/5 + 1/5))/cosh(Cx1)^2 - ((tanh(Cx1)^2 - ((9*exp(-(2*t)/5))/10 + 1/10)^2)*((9*exp(-(2*t)/5))/10 + 1/10))/cosh(Cx1)^2)*(tanh(Cx1)^2 - ((9*exp(-(2*t)/5))/10 + 1/10)^2)^2) - (tanh(Cx1)*((2*exp(-(2*t)/5)*tanh(Cx1))/5 - (4*exp(-(2*t)/5)*tanh(Cx1)*((9*exp(-(2*t)/5))/10 + 1/10)^2)/5)*((9*exp(-(2*t)/5))/10 + 1/10)*((2*tanh(Cx1)*(tanh(Cx1)^2 - 1)*((9*exp(-(2*t)/5))/5 + 1/5))/cosh(Cx1)^2 - (2*tanh(Cx1)*(tanh(Cx1)^2 - 1)*((9*exp(-(2*t)/5))/10 + 1/10))/cosh(Cx1)^2 + (2*sinh(Cx1)*tanh(Cx1)^2*((9*exp(-(2*t)/5))/5 + 1/5))/cosh(Cx1)^3 - (2*sinh(Cx1)*(tanh(Cx1)^2 - ((9*exp(-(2*t)/5))/10 + 1/10)^2)*((9*exp(-(2*t)/5))/10 + 1/10))/cosh(Cx1)^3))/(((tanh(Cx1)^2*((9*exp(-(2*t)/5))/5 + 1/5))/cosh(Cx1)^2 - ((tanh(Cx1)^2 - ((9*exp(-(2*t)/5))/10 + 1/10)^2)*((9*exp(-(2*t)/5))/10 + 1/10))/cosh(Cx1)^2)^2*(tanh(Cx1)^2 - ((9*exp(-(2*t)/5))/10 + 1/10)^2)) + (Cthetahat*Cx1*(tanh(Cx1)^2 - ((9*exp(-(2*t)/5))/10 + 1/10)^2)*(tanh(Cx1)^2 - 1))/(tanh(Cx1)^2*((9*exp(-(2*t)/5))/10 + 1/10))))/(tanh(Cx1)^2 - ((9*exp(-(2*t)/5))/10 + 1/10)^2);    
    
    Cdalpha1_dbeta = (tanh(Cx1)*((tanh(Cx1)^2 - Cbeta^2)^2/((2*Cbeta*tanh(Cx1)^2)/cosh(Cx1)^2 - (Cbeta*(tanh(Cx1)^2 - Cbeta^2))/cosh(Cx1)^2) + (k1*(tanh(Cx1)^2 - Cbeta^2)^2)/((2*Cbeta*tanh(Cx1)^2)/cosh(Cx1)^2 - (Cbeta*(tanh(Cx1)^2 - Cbeta^2))/cosh(Cx1)^2) - (Cx1*(tanh(Cx1)^2 - Cbeta^2))/(Cbeta*tanh(Cx1)) - (Cthetahat*Cx1*(tanh(Cx1)^2 - Cbeta^2))/(Cbeta*tanh(Cx1)) + (Cbeta*tanh(Cx1)*((2*exp(-(2*t)/5)*tanh(Cx1))/5 - (4*Cbeta^2*exp(-(2*t)/5)*tanh(Cx1))/5))/(((2*Cbeta*tanh(Cx1)^2)/cosh(Cx1)^2 - (Cbeta*(tanh(Cx1)^2 - Cbeta^2))/cosh(Cx1)^2)*(tanh(Cx1)^2 - Cbeta^2))))/(tanh(Cx1)^2 - Cbeta^2) + (Cbeta*tanh(Cx1)*((2*Cx1)/tanh(Cx1) - (4*Cbeta*(tanh(Cx1)^2 - Cbeta^2))/((2*Cbeta*tanh(Cx1)^2)/cosh(Cx1)^2 - (Cbeta*(tanh(Cx1)^2 - Cbeta^2))/cosh(Cx1)^2) - ((tanh(Cx1)^2 - Cbeta^2)^2*((2*Cbeta^2)/cosh(Cx1)^2 - (tanh(Cx1)^2 - Cbeta^2)/cosh(Cx1)^2 + (2*tanh(Cx1)^2)/cosh(Cx1)^2))/((2*Cbeta*tanh(Cx1)^2)/cosh(Cx1)^2 - (Cbeta*(tanh(Cx1)^2 - Cbeta^2))/cosh(Cx1)^2)^2 + (2*Cthetahat*Cx1)/tanh(Cx1) + (tanh(Cx1)*((2*exp(-(2*t)/5)*tanh(Cx1))/5 - (4*Cbeta^2*exp(-(2*t)/5)*tanh(Cx1))/5))/(((2*Cbeta*tanh(Cx1)^2)/cosh(Cx1)^2 - (Cbeta*(tanh(Cx1)^2 - Cbeta^2))/cosh(Cx1)^2)*(tanh(Cx1)^2 - Cbeta^2)) + (Cx1*(tanh(Cx1)^2 - Cbeta^2))/(Cbeta^2*tanh(Cx1)) - (4*Cbeta*k1*(tanh(Cx1)^2 - Cbeta^2))/((2*Cbeta*tanh(Cx1)^2)/cosh(Cx1)^2 - (Cbeta*(tanh(Cx1)^2 - Cbeta^2))/cosh(Cx1)^2) - (k1*(tanh(Cx1)^2 - Cbeta^2)^2*((2*Cbeta^2)/cosh(Cx1)^2 - (tanh(Cx1)^2 - Cbeta^2)/cosh(Cx1)^2 + (2*tanh(Cx1)^2)/cosh(Cx1)^2))/((2*Cbeta*tanh(Cx1)^2)/cosh(Cx1)^2 - (Cbeta*(tanh(Cx1)^2 - Cbeta^2))/cosh(Cx1)^2)^2 + (Cthetahat*Cx1*(tanh(Cx1)^2 - Cbeta^2))/(Cbeta^2*tanh(Cx1)) + (2*Cbeta^2*tanh(Cx1)*((2*exp(-(2*t)/5)*tanh(Cx1))/5 - (4*Cbeta^2*exp(-(2*t)/5)*tanh(Cx1))/5))/(((2*Cbeta*tanh(Cx1)^2)/cosh(Cx1)^2 - (Cbeta*(tanh(Cx1)^2 - Cbeta^2))/cosh(Cx1)^2)*(tanh(Cx1)^2 - Cbeta^2)^2) - (8*Cbeta^2*exp(-(2*t)/5)*tanh(Cx1)^2)/(5*((2*Cbeta*tanh(Cx1)^2)/cosh(Cx1)^2 - (Cbeta*(tanh(Cx1)^2 - Cbeta^2))/cosh(Cx1)^2)*(tanh(Cx1)^2 - Cbeta^2)) - (Cbeta*tanh(Cx1)*((2*exp(-(2*t)/5)*tanh(Cx1))/5 - (4*Cbeta^2*exp(-(2*t)/5)*tanh(Cx1))/5)*((2*Cbeta^2)/cosh(Cx1)^2 - (tanh(Cx1)^2 - Cbeta^2)/cosh(Cx1)^2 + (2*tanh(Cx1)^2)/cosh(Cx1)^2))/(((2*Cbeta*tanh(Cx1)^2)/cosh(Cx1)^2 - (Cbeta*(tanh(Cx1)^2 - Cbeta^2))/cosh(Cx1)^2)^2*(tanh(Cx1)^2 - Cbeta^2))))/(tanh(Cx1)^2 - Cbeta^2) + (2*Cbeta^2*tanh(Cx1)*((tanh(Cx1)^2 - Cbeta^2)^2/((2*Cbeta*tanh(Cx1)^2)/cosh(Cx1)^2 - (Cbeta*(tanh(Cx1)^2 - Cbeta^2))/cosh(Cx1)^2) + (k1*(tanh(Cx1)^2 - Cbeta^2)^2)/((2*Cbeta*tanh(Cx1)^2)/cosh(Cx1)^2 - (Cbeta*(tanh(Cx1)^2 - Cbeta^2))/cosh(Cx1)^2) - (Cx1*(tanh(Cx1)^2 - Cbeta^2))/(Cbeta*tanh(Cx1)) - (Cthetahat*Cx1*(tanh(Cx1)^2 - Cbeta^2))/(Cbeta*tanh(Cx1)) + (Cbeta*tanh(Cx1)*((2*exp(-(2*t)/5)*tanh(Cx1))/5 - (4*Cbeta^2*exp(-(2*t)/5)*tanh(Cx1))/5))/(((2*Cbeta*tanh(Cx1)^2)/cosh(Cx1)^2 - (Cbeta*(tanh(Cx1)^2 - Cbeta^2))/cosh(Cx1)^2)*(tanh(Cx1)^2 - Cbeta^2))))/(tanh(Cx1)^2 - Cbeta^2)^2;
    Cdalpha1_dthetahat = -Cx1; 
    CW2 = [-Cdalpha1_dx1*CW1;
                         0];
    Cpsibar = [CPi+  Cdalpha1_dx1*Ckappa1 - Cdalpha1_dbeta * Cdbeta /Cz1   - Cdalpha1_dthetahat*gamma*CPi*Cx1 - Cdalpha1_dx1*CW1*Cthetahat; 
              -Cdalpha1_dx1+Cdalpha1_dthetahat*gamma*Cdalpha1_dx1*Cx1];       
    Cz2 = Cx2 - Calpha1;
    Ckappa2 = k2+0.5*( delta * (norm(CW2,2))  + delta + 1 + (norm(Cpsibar,2))  );  
    Cubar = -Ckappa2 * Cz2;
    Cu = Crhohat * Cubar; 
    Cdthetahat = gamma * CPi * Cx1 * Cz1  - gamma * Cdalpha1_dx1 * Cx1 *Cz2;
    Cthetahat = Cthetahat + Cdthetahat * step; 
    Cdrhohat = -gamma  * Cz2 * Cubar;
    Crhohat = Crhohat + Cdrhohat * step;  
    Cdx1 = Cx2+Cthetat*Cx1;
    Cdx2 = Cbt*Cu;
    Cx1 = Cx1 + Cdx1 * step;    
    Cx2 = Cx2 + Cdx2 * step;   
    Cx1s(i)=Cx1;
    Cx2s(i)=Cx2;
    Cus(i) = Cu;
    Cbs(i)=Cboundary;
    Cbs1(i)=-Cboundary; 
    Cthetas(i)=Cthetahat;
    Crhos(i)=Crhohat;
    Cthetats(i)=Cthetat;
    Cbts(i)=Cbt; 
end

%%   画图------------------------------chen--------------ye-----------------chen---------------ye-----------
figure(1)   
plot(ts,Ax1s,ts,Bx1s,ts,Cx1s, ts,Bbs,ts,Bbs1,  ts,Cbs,ts,Cbs1, 'linewidth',1);  
legend({'$Controller~1^{[1]}$','$Controller~2$','$Controller~3$'},'interpreter', 'latex');  
ylabel('$y(t)$','interpreter','latex','linewidth',18);
xlabel('Time(s)','interpreter','latex','linewidth',18);
grid on;  

figure(2)   
plot(ts,Ax2s,ts,Bx2s,ts,Cx2s,  'linewidth',1);  
legend({'$Controller~1^{[1]}$','$Controller~2$','$Controller~3$'},'interpreter', 'latex');  
ylabel('$x_2(t)$','interpreter','latex','linewidth',18);
xlabel('Time(s)','interpreter','latex','linewidth',18);
grid on;  

figure(3) 
plot(ts,Aus,ts,Bus,ts,Cus,  'linewidth',1);  
legend({'$Controller~1^{[1]}$','$Controller~2$','$Controller~3$'},'interpreter', 'latex');  
xlabel('Time(s)','interpreter','latex','linewidth',18);
ylabel('Control Input $u(t)$','interpreter','latex','linewidth',18);
grid on; 
axes('Position',[0.4 0.4 0.45 0.15],'XGrid','on','XTickLabel',[])        %生成子图  0.4距左边距离   0.2距下边距离    0.4是小图长度  0.2是小图宽度 
plot(ts,Aus,ts,Bus,ts,Cus, 'linewidth',1);
axis([0 3 -40 30]); 


figure(4)
plot(ts,Athetas,ts,Bthetas,ts,Cthetas,  'linewidth',1);  
legend({'$Controller~1^{[1]}$','$Controller~2$','$Controller~3$'},'interpreter', 'latex');  
xlabel('Time(s)','interpreter','latex','linewidth',18);
ylabel('$\hat{\theta}(t) $','interpreter','latex','linewidth',18);
grid on; 

figure(5) 
plot(ts,Arhos,ts,Brhos,ts,Crhos,  'linewidth',1);  
legend({'$Controller~1^{[1]}$','$Controller~2$','$Controller~3$'},'interpreter', 'latex');   
xlabel('Time(s)','interpreter','latex','linewidth',18);
ylabel('$\hat{\rho}(t) $','interpreter','latex','linewidth',18);
grid on;

figure(6)
subplot(2,1,1)
plot(ts,Athetats,ts,Bthetats,ts,Cthetats,  'linewidth',1);  
legend({'$Controller~1^{[1]}$','$Controller~2$','$Controller~3$'},'interpreter', 'latex');  
xlabel('Time(s)','interpreter','latex','linewidth',18);
ylabel('$\theta(t)$ ','interpreter','latex','linewidth',18);
grid on; 
 
subplot(2,1,2)
plot(ts,Bbts,ts,Cbts,ts,Abts,  'linewidth',1);  
legend({'$Controller~2$','$Controller~3$','$Controller~1^{[1]}$'},'interpreter', 'latex'); 
xlabel('Time(s)','interpreter','latex','linewidth',18);
ylabel('$b(t)$ ','interpreter','latex','linewidth',18);
grid on;  
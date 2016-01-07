% Written by mhamad mahdi salati

% investigation for axial loads
%    The Plane Frame Element K M  Global Matrises
% static anlysis fr step load ,impact load and cretical load
% dynamic anlysis for impact load and harmnic load and...
% euler bernoli beam for mass matrix and stiffness matrix and damping matrix
% newton raphson metod for analysis step gooz
% newmark family metod for analysis dynamic
% impact : impacted mass by bar
%               : impact load
% all by fineit element metod


function []=static_and_dynamic()
clc
% [Node ncon Constraints]= input6();
global A E Ks G I sncon ro_t nsec dx ncon Node sNode_q ax ay az am dt Ty
flag = menu('Select',...
						'Dynamic',...
						'static',...
                        'static of dynamic');
                    if  flag==2;
                         flag1 = 2;%non damping for analysis static
                    elseif flag==1
flag1 = menu('Select',...
						'daming',...
						'no daming');
disp('Please wait ...');
                    else
                        flag1 =2;%non damping for analysis static
                    end
%% Preliminaries of the main while loop
% Deta 0
% A=1; E=30e6;  G=12e6;   I=1/12;
% Deta1
%% deta 1
E= 29e9;    ro_t= 2500;     nu=0.2;     bbbb=0.5;   hhhh=0.3;   L=6;
M2=ro_t*bbbb*hhhh*L;    M1=45*M2;
v0=0;   max_iter=100;
%% deta 2
% E= 210e9;     ro_t=7800;      stress_cr=1500e6;       L=0.193;    nu=0.2;       
% bbbb=8.7*10^-3;     hhhh=0.63*10^-3;    M2=ro_t*bbbb*hhhh*L;    M1=0.369;
% v0=0;   max_iter=50;
%% calc
nsec=10;%Number elemenets beam
Node=zeros(nsec+1,3);
Node(:,1)=0:L/nsec:L;
ncon=zeros(nsec,7);
ncon(:,1)=1:nsec;
ncon(:,2)=2:nsec+1;
c1=1e27;
c2=1e-27;
Constraints=[1	c2	c1	c1
    nsec+1	c1	c1	c1];

A=bbbb*hhhh;        I=bbbb*hhhh^3/12;       r=sqrt(I/A);        ct=sqrt(E/ro_t);
nx=nsec+1;      dx=L/nsec;      %dt=dx/ct;       
% dt=dx/10/ct;

Ty=L/ct;        
TT=5*Ty;
dt=Ty/8;
if flag1 == 1
coef_damp=0.00003;
else
coef_damp=0;%damping
end
nt=floor(TT/dt)+1;% step time
G=E/(2*(1+nu));
%% stiffness Deta
Ks=5/6;
% Ks=1; % Ks is applied in ay and az.
Axx=E*A;
Dxx=E*I;
Sxx=Ks*G*A*0;

%% Deta matrix
ax=1;       ay=1;       az=1;       
if flag == 1
    am=1;
else
    am=0;
end

sncon=size(ncon,1);
sNode=size(Node,1);   sNode_q=2*sNode-1;
wabc=zeros(sncon,3);      uabc=zeros(sncon,3);
% K_sys=zeros(3*sNode,3*sNode);    F_sys=zeros(3*sNode,1);

X_sys_former=zeros(3*sNode_q,1);
x_ele=zeros(sncon,9);  %10*ones(sncon,6);
%%  Applying Boundary Conditions
c1=1e270;
c2=1e-270;
%[node_number rigidity_u rigidity_w rigidity_teta] 
% Constraints=[1	c2	c1	c2
% 5	c1	c2	c1];
Constraints=[1	c2	c1	c2
nsec+1	c1	c1	c2];
[as bs]=find(Constraints>1000);
np_deg=sort(2*3*(-1+Constraints(as))+bs-ones(size(bs,1),1));
% ----------------------
% wab_former=zeros(sncon,1);
%% LOADING static and Deta
f0=0;
q=0.0001;  
load_increment_q0=0.0001;

F0_=3500;  
LL=100;

load_increment_F0=350;

NI_F0=abs(F0_/load_increment_F0);
inp_graph=zeros(NI_F0,1);
F0=load_increment_F0;
Determinant=zeros(NI_F0,1);

NI_q0=q/load_increment_q0;
q0=load_increment_q0;
X_load_incre=zeros(3*sNode_q,NI_F0);

F_app=0;
if flag == 1
dynamic(nt,v0,coef_damp,np_deg,max_iter,F0,x_ele,f0,q0);
elseif flag == 2
 static(F_app,x_ele,f0,q0,X_load_incre,NI_F0,np_deg,load_increment_F0);
else
  dynamic(nt,v0,coef_damp,np_deg,max_iter,F0,x_ele,f0,q0);
end

 %% function uw_element
function [wabc,uabc,x_ele]=uw_ele(ncon,X_sys,sncon)

for el=1:sncon
  n1=ncon(el,1);
  n2=ncon(el,2);
%   r1=3*(n1-1)+1:3*(n1-1)+3;
%   r2=3*(n2-1)+1:3*(n2-1)+3;
  
  r1=3*(n1-1+el-1)+1:3*(n1-1+el-1)+3;   % Note: Just for Unbranched frams
  r2=3*(n2-1+el-1)+1:3*(n2-1+el-1)+3;
  r3=3*(n2-1+el)+1:3*(n2-1+el)+3;
 
  rr=[r1 r2 r3];  
  x_ele(el,:)=X_sys(rr);

  wabc(el,:)=x_ele(el,[2 5 8]);
  
  uabc(el,:)=x_ele(el,[1 4 7]);
end

%% function calculate_K_T_element
function [K_ele, T_ele]=calculate_K_T_element(xa,xc,A,E,Ks,G,I,wabc_e,uabc_e)
Axx=E*A;
Dxx=E*I;
Sxx=Ks*G*A;
xb=0.5*(xa+xc);

i=1:3 ; j=1:3;

NGP=3;
LGP=2;
Gauss_full=gauss(NGP);
Gauss_reduced=gauss(LGP);

%functions needed to be implemented by full integration
% f_K11=(DSx(i))'*DSx(j);
F_K11=0;
for ite=1:NGP
    xi=Gauss_full(ite,1);  wt=Gauss_full(ite,2);
%     S1 =[ 1/2-1/2*xi, 1/2+1/2*xi];
   Jac=xa*xi-1/2*xa-2*xb*xi+xc*xi+1/2*xc;   % Jac =-1/2*xa+1/2*xb;
   DSx =[ -(2*xi-1)/(-2*xa*xi+xa+4*xb*xi-2*xc*xi-xc),      4*xi/(-2*xa*xi+xa+4*xb*xi-2*xc*xi-xc), -(2*xi+1)/(-2*xa*xi+xa+4*xb*xi-2*xc*xi-xc)];
 %DSx =[1/(xa-xb), -1/(xa-xb)];
    
    f_K11=(DSx(i))'*DSx(j);
    F_K11=F_K11+f_K11*wt*Jac;
end
 K11=Axx*F_K11;
 %functions needed to be implemented by reduced integration

wa=wabc_e(1);    wb=wabc_e(2);    wc=wabc_e(3);
ua=uabc_e(1);      ub=uabc_e(2);     uc=uabc_e(3);
F_K12=0;
F_K22_1=0;
F_K22_2=0;
F_K23=0;
F_K33_1=0;
F_K33_2=0;
F_T12_extra=0;    % the word extra is introduced here because T12 has other terms: K12
F_T22_extra_1=0;
F_T22_extra_2=0;
for ite=1:LGP
    xi=Gauss_reduced(ite,1);  wt=Gauss_reduced(ite,2);
    S2(1)=-0.5*(1-xi)*xi;  S2(2)=(1-xi*xi);    S2(3)=0.5*(1+xi)*xi;   %S1 =[ 1/2-1/2*xi, 1/2+1/2*xi];
    Jac=xa*xi-1/2*xa-2*xb*xi+xc*xi+1/2*xc;   % Jac =-1/2*xa+1/2*xb;
   DSx =[ -(2*xi-1)/(-2*xa*xi+xa+4*xb*xi-2*xc*xi-xc),      4*xi/(-2*xa*xi+xa+4*xb*xi-2*xc*xi-xc), -(2*xi+1)/(-2*xa*xi+xa+4*xb*xi-2*xc*xi-xc)];
 %DSx =[1/(xa-xb), -1/(xa-xb)];
%     f_T12_extra=(wa*DSx(1)+wb*DSx(2)).*((DSx(i))'*DSx(j));
%     F_T12_extra=F_T12_extra+f_T12_extra*wt*Jac;
    
    f_T22_extra_1=(ua*DSx(1)+ub*DSx(2)+uc*DSx(3)).*((DSx(i))'*DSx(j));
    F_T22_extra_1=F_T22_extra_1+f_T22_extra_1*wt*Jac;
    
    f_T22_extra_2=((wa*DSx(1)+wb*DSx(2)+wc*DSx(3))^2).*(DSx(i))'*DSx(j);
    F_T22_extra_2=F_T22_extra_2+f_T22_extra_2*wt*Jac;
    
    f_K12=(wa*DSx(1)+wb*DSx(2)+wc*DSx(3)).*((DSx(i))'*DSx(j));
    F_K12=F_K12+f_K12*wt*Jac;
    
    f_K22_1=(DSx(i))'*DSx(j);
    F_K22_1=F_K22_1+f_K22_1*wt*Jac;
    
    f_K22_2=((wa*DSx(1)+wb*DSx(2)+wc*DSx(3))^2).*(DSx(i))'*DSx(j);
    F_K22_2=F_K22_2+f_K22_2*wt*Jac;
    
    f_K23=(DSx(i))'*S2(j);   % or maybe this is true:   DSx(i)*(S1(j))';
    F_K23=F_K23+f_K23*wt*Jac;
    
    f_K33_1=(DSx(i))'*DSx(j);
    F_K33_1=F_K33_1+f_K33_1*wt*Jac;
    
    f_K33_2=(S2(i))'*S2(j);
    F_K33_2=F_K33_2+f_K33_2*wt*Jac;    
end
%  K11=Axx*F_K11;
%  K12=0.5*Axx*F_K12;   K21=2*K12;
 K12=0;   K21=2*K12;
 K13=zeros(3,3);    K31=K13;
 K22=Sxx*F_K22_1+0.5*Axx*F_K22_2;
 K23=Sxx*F_K23;     K32=K23';
 K33=Dxx*F_K33_1+Sxx*F_K33_2;

 K_ele([1 4 7],[1 4 7])=K11;
 K_ele([2 5 8],[2 5 8])=K22;
 K_ele([3 6 9],[3 6 9])=K33;
 K_ele([2 5 8],[1 4 7])=K21;
 K_ele([1 4 7],[2 5 8])=K12;
 K_ele([1 4 7],[3 6 9])=K13;
 K_ele([2 5 8],[3 6 9])=K23;
 K_ele([3 6 9],[1 4 7])=K31;
 K_ele([3 6 9],[2 5 8])=K32;
 
 T11=K11;
 T12=K12+0.5*Axx*F_T12_extra;    %T12=2*K12;
%  T12=2*K12;
 T21=K21;
 T13=K13;    
 T31=K31;
 T22=K22+Axx*(F_T22_extra_1+F_T22_extra_2);
 T23=K23;     
 T32=K32;
 T33=K33;
 
 T_ele([1 4 7],[1 4 7])=T11;
 T_ele([2 5 8],[2 5 8])=T22;
 T_ele([3 6 9],[3 6 9])=T33;
 T_ele([2 5 8],[1 4 7])=T21;
 T_ele([1 4 7],[2 5 8])=T12;
 T_ele([1 4 7],[3 6 9])=T13;
 T_ele([2 5 8],[3 6 9])=T23;
 T_ele([3 6 9],[1 4 7])=T31;
 T_ele([3 6 9],[2 5 8])=T32;
dfgf=0;
return
 
%% function calculate_K__element 3NODE
function [K_ele]=calculate_K__element(xa,xc,A,E,Ks,G,I,wabc_e,uabc_e,ax,ay,az)
Axx=E*A;
Dxx=E*I;
Sxx=Ks*G*A;
xb=0.5*(xa+xc);

i=1:3 ; j=1:3;

NGP=2;
LGP=2;
Gauss_full=gauss(NGP);
Gauss_reduced=gauss(LGP);

%functions needed to be implemented by full integration
% f_K11=(DSx(i))'*DSx(j);
F_K11=0;
for ite=1:NGP
    xi=Gauss_full(ite,1);  wt=Gauss_full(ite,2);
    %     S1 =[ 1/2-1/2*xi, 1/2+1/2*xi];
    Jac=xa*xi-1/2*xa-2*xb*xi+xc*xi+1/2*xc;   % Jac =-1/2*xa+1/2*xb;
    DSx =[ -(2*xi-1)/(-2*xa*xi+xa+4*xb*xi-2*xc*xi-xc),      4*xi/(-2*xa*xi+xa+4*xb*xi-2*xc*xi-xc), -(2*xi+1)/(-2*xa*xi+xa+4*xb*xi-2*xc*xi-xc)];
    %DSx =[1/(xa-xb), -1/(xa-xb)];
    
    f_K11=ax*(DSx(i))'*DSx(j);
    F_K11=F_K11+f_K11*wt*Jac;
end
K11=Axx*F_K11;
%functions needed to be implemented by reduced integration

wa=wabc_e(1);    wb=wabc_e(2);    wc=wabc_e(3);
ua=uabc_e(1);      ub=uabc_e(2);     uc=uabc_e(3);
F_K12=0;
F_K22_1=0;
F_K22_2=0;
F_K23=0;
F_K33_1=0;
F_K33_2=0;

for ite=1:LGP
    xi=Gauss_reduced(ite,1);  wt=Gauss_reduced(ite,2);
    S2(1)=-0.5*(1-xi)*xi;  S2(2)=(1-xi*xi);    S2(3)=0.5*(1+xi)*xi;   %S1 =[ 1/2-1/2*xi, 1/2+1/2*xi];
    Jac=xa*xi-1/2*xa-2*xb*xi+xc*xi+1/2*xc;   % Jac =-1/2*xa+1/2*xb;
    DSx =[ -(2*xi-1)/(-2*xa*xi+xa+4*xb*xi-2*xc*xi-xc),      4*xi/(-2*xa*xi+xa+4*xb*xi-2*xc*xi-xc), -(2*xi+1)/(-2*xa*xi+xa+4*xb*xi-2*xc*xi-xc)];
    
    f_K12=ax*ay*(wa*DSx(1)+wb*DSx(2)+wc*DSx(3)).*((DSx(i))'*DSx(j));
    F_K12=F_K12+f_K12*wt*Jac;
    
    f_K22_1=az*(DSx(i))'*DSx(j);
    F_K22_1=F_K22_1+f_K22_1*wt*Jac;
    
    f_K22_2=ax*ay*((wa*DSx(1)+wb*DSx(2)+wc*DSx(3))^2).*(DSx(i))'*DSx(j);
    F_K22_2=F_K22_2+f_K22_2*wt*Jac;
    
    f_K23=az*(DSx(i))'*S2(j);   % or maybe this is true:   DSx(i)*(S1(j))';
    F_K23=F_K23+f_K23*wt*Jac;
    
    f_K33_1=ay*(DSx(i))'*DSx(j);
    F_K33_1=F_K33_1+f_K33_1*wt*Jac;
    
    f_K33_2=az*(S2(i))'*S2(j);
    F_K33_2=F_K33_2+f_K33_2*wt*Jac;
end
%  K11=Axx*F_K11;
K12=0.5*Axx*F_K12;   K21=2*K12;
K13=zeros(3,3);    K31=K13;
K22=Sxx*F_K22_1+0.5*Axx*F_K22_2;
K23=Sxx*F_K23;     K32=K23';
K33=Dxx*F_K33_1+Sxx*F_K33_2;

K_ele([1 4 7],[1 4 7])=K11;
K_ele([2 5 8],[2 5 8])=K22;
K_ele([3 6 9],[3 6 9])=K33;
K_ele([2 5 8],[1 4 7])=K21;
K_ele([1 4 7],[2 5 8])=K12;
K_ele([1 4 7],[3 6 9])=K13;
K_ele([2 5 8],[3 6 9])=K23;
K_ele([3 6 9],[1 4 7])=K31;
K_ele([3 6 9],[2 5 8])=K32;

return

%% function calculate_K_element 2NODE
function K_ele=calculate_K_element(xa,xb,A,E,Ks,G,I,wab_e)
Axx=E*A;
Dxx=E*I;
Sxx=Ks*G*A;
xab=[xa xb];

i=1:2 ; j=1:2;

NGP=2;
LGP=1;
Gauss_full=gauss(NGP);
Gauss_reduced=gauss(LGP);

%functions needed to be implemented by full integration
% f_K11=(DSx(i))'*DSx(j);
F_K11=0;
for ite=1:NGP
    xi=Gauss_full(ite,1);  wt=Gauss_full(ite,2);
%     S1 =[ 1/2-1/2*xi, 1/2+1/2*xi];
    Jac =-1/2*xa+1/2*xb;
    DSx =[1/(xa-xb), -1/(xa-xb)];
    
    f_K11=(DSx(i))'*DSx(j);
    F_K11=F_K11+f_K11*wt*Jac;
end
 K11=Axx*F_K11;
 %functions needed to be implemented by reduced integration

wa=wab_e(1);    wb=wab_e(2);
F_K12=0;
F_K22_1=0;
F_K22_2=0;
F_K23=0;
F_K33_1=0;
F_K33_2=0;
for ite=1:LGP
    xi=Gauss_reduced(ite,1);  wt=Gauss_reduced(ite,2);
    S1 =[ 1/2-1/2*xi, 1/2+1/2*xi];
    Jac =-1/2*xa+1/2*xb;
    DSx =[1/(xa-xb), -1/(xa-xb)];
    
    f_K12=(wa*DSx(1)+wb*DSx(2)).*((DSx(i))'*DSx(j));
    F_K12=F_K12+f_K12*wt*Jac;
    
    f_K22_1=(DSx(i))'*DSx(j);
    F_K22_1=F_K22_1+f_K22_1*wt*Jac;
    
    f_K22_2=((wa*DSx(1)+wb*DSx(2))^2).*(DSx(i))'*DSx(j);
    F_K22_2=F_K22_2+f_K22_2*wt*Jac;
    
    f_K23=(DSx(i))'*S1(j);   % or maybe this is true:   DSx(i)*(S1(j))';
    F_K23=F_K23+f_K23*wt*Jac;
    
    f_K33_1=(DSx(i))'*DSx(j);
    F_K33_1=F_K33_1+f_K33_1*wt*Jac;
    
    f_K33_2=(S1(i))'*S1(j);
    F_K33_2=F_K33_2+f_K33_2*wt*Jac;    
end
 K11=Axx*F_K11;
 K12=0.5*Axx*F_K12;   K21=2*K12;
 K13=zeros(2,2);    K31=K13;
 K22=Sxx*F_K22_1+0.5*Axx*F_K22_2;
 K23=Sxx*F_K23;     K32=K23';
 K33=Dxx*F_K33_1+Sxx*F_K33_2;
%  
%  K32(1,2)=- K32(1,2);
%   K32(2,1)=- K32(2,1);
%  
 
%  ELK=zeros(6,6);
 ELK([1 4],[1 4])=K11;
 ELK([2 5],[2 5])=K22;
 ELK([3 6],[3 6])=K33;
 ELK([2 5],[1 4])=K21;
 ELK([1 4],[2 5])=K12;
 ELK([1 4],[3 6])=K13;
 ELK([2 5],[3 6])=K23;
 ELK([3 6],[1 4])=K31;
 ELK([3 6],[2 5])=K32;

 K_ele=ELK;    % K_ele=double(ELK);
 return
 
%% function Mass
function [M_ele]=calculate_M_element(xa,xc,ro_t,A,I)

RoA=ro_t*A;
RoI=ro_t*A;
% RoI=ro_t*I;
xb=0.5*(xa+xc);

i=1:3 ; j=1:3;

NGP=3;
LGP=2;
Gauss_full=gauss(NGP);
Gauss_reduced=gauss(LGP);

%functions needed to be implemented by full integration
% f_K11=(DSx(i))'*DSx(j);
F_M11=0;
for ite=1:NGP
    xi=Gauss_full(ite,1);  wt=Gauss_full(ite,2);
    
    S2(1)=-0.5*(1-xi)*xi;  S2(2)=(1-xi*xi);    S2(3)=0.5*(1+xi)*xi;   %S1 =[ 1/2-1/2*xi, 1/2+1/2*xi];
    
    Jac=xa*xi-1/2*xa-2*xb*xi+xc*xi+1/2*xc;   % Jac =-1/2*xa+1/2*xb;
    
    f_M11=(S2(i))'*S2(j);%(DSx(i))'*DSx(j);
    F_M11=F_M11+f_M11*wt*Jac;
    
end
M11=RoA*F_M11;
%functions needed to be implemented by reduced integration

F_M33=0;

for ite=1:LGP
    xi=Gauss_reduced(ite,1);  wt=Gauss_reduced(ite,2);
    S2(1)=-0.5*(1-xi)*xi;  S2(2)=(1-xi*xi);    S2(3)=0.5*(1+xi)*xi;   %S1 =[ 1/2-1/2*xi, 1/2+1/2*xi];
    Jac=xa*xi-1/2*xa-2*xb*xi+xc*xi+1/2*xc;   % Jac =-1/2*xa+1/2*xb;
    f_M33=(S2(i))'*S2(j);
    F_M33=F_M33+f_M33*wt*Jac;
end
F_M22=F_M33;;
M22=RoA*F_M22;
M33=RoI*F_M33;


M13=zeros(3,3);    M31=M13;
M23=zeros(3,3);    M32=M23;
M12=zeros(3,3);    M21=M12;

M_ele([1 4 7],[1 4 7])=M11;
M_ele([2 5 8],[2 5 8])=M22;
M_ele([3 6 9],[3 6 9])=M33;
M_ele([2 5 8],[1 4 7])=M21;
M_ele([1 4 7],[2 5 8])=M12;
M_ele([1 4 7],[3 6 9])=M13;
M_ele([2 5 8],[3 6 9])=M23;
M_ele([3 6 9],[1 4 7])=M31;
M_ele([3 6 9],[2 5 8])=M32;
return

%% function Mass
function [m]=mass_mms(xa,xc,ro,A,I,dx);
    
xb=xc-xa;
m11=0;
for i=1:3
    for j=1:3
    syms x
%     S1 =[ 1/2-1/2*xi, 1/2+1/2*xi];
    S2(1)=(2*x^2)-(3*x)+1;  S2(2)=(-4*x^2)+(4*x);    S2(3)=(2*x^2)-x;   %S1 =[ 1-x, x];
    m11(i,j)=int(S2(i)*S2(j),0,1);
    end
end
 %functions needed to be implemented by reduced integration
 
 S3=[(x^2)-(5*(x^3)/4)-((x^4)/2)+(3*(x^5))
    (dx/4)*(x^2-x^3-x^4+x^5)
    1-(2*x^2)+x^4
    dx*(x-(2*x^2)+x^5)
    x^2+(5*(x^3)/4)-((x^4)/2)-(3*(x^5)/4)
    (dx/4)*(-x^2-x^3+x^4+x^5)];

    syms x
for i=1:6
    for j=1:6
    m22(i,j)=int(S3(i)*S3(j),0,1);
    end
end
 m11=m11*ro*A;
 m22=m22;
 m_ele([1 4 7],[1 4 7])=m11;
 m_ele([2 5 8],[2 5 8])=m22([1 3 5],[1 3 5])*ro*A;
 m_ele([3 6 9],[3 6 9])=m22([2 4 6],[2 4 6])*ro*I;
 m_ele([2 5 8],[1 4 7])=0;%m21;
 m_ele([1 4 7],[2 5 8])=0;%m12;
 m_ele([1 4 7],[3 6 9])=0;%m13;
 m_ele([2 5 8],[3 6 9])=0;%m22([1 3 5],[2 4 6]);
 m_ele([3 6 9],[1 4 7])=0;%m31;
 m_ele([3 6 9],[2 5 8])=0;%m22([2 4 6],[1 3 5]);
 m= m_ele;
 hhh=54;

%% function  calculate_F_element
function F_ele=calculate_F_element(xa,xc,f0,q0,F0)

xb=0.5*(xa+xc);
i=1:3 ; j=1:3;

NGP=3;
LGP=2;
Gauss_full=gauss(NGP);
Gauss_reduced=gauss(LGP);

%functions needed to be implemented by full integration
F_F1=0;
F_F2=0;
for ite=1:NGP
    xi=Gauss_full(ite,1);  wt=Gauss_full(ite,2);
    S2(1)=-0.5*(1-xi)*xi;  S2(2)=(1-xi*xi);    S2(3)=0.5*(1+xi)*xi;   %S1 =[ 1/2-1/2*xi, 1/2+1/2*xi];
    Jac=xa*xi-1/2*xa-2*xb*xi+xc*xi+1/2*xc;   % Jac =-1/2*xa+1/2*xb;
%    DSx =[ -(2*xi-1)/(-2*xa*xi+xa+4*xb*xi-2*xc*xi-xc),      4*xi/(-2*xa*xi+xa+4*xb*xi-2*xc*xi-xc), -(2*xi+1)/(-2*xa*xi+xa+4*xb*xi-2*xc*xi-xc)];
%     DSx =[1/(xa-xb), -1/(xa-xb)];
    f_F1=S2(i);
    F_F1=F_F1+f_F1*wt*Jac;
    f_F2=S2(i);        %S2(i);
    F_F2=F_F2+f_F2*wt*Jac;
end
%  K11=Axx*F_K11;
 F1=f0*F_F1;
 F2=q0*F_F2;
 F3=zeros(3,1);    %F3=q0*F_F3;
F_ele([1;4;7])=F1;
F_ele([2;5;8])=F2;
F_ele([3;6;9])=F3;
F_ele=F_ele';     %F_ele=double(F_ele');
% F_ele=vpa(F_ele');
if xa==0
    F_ele(1,1)=F0;
end
 return

 %% function  Assembling_K_F
function [K_sys,F_sys]=K_F_sys(sNode_q,Node,ncon,uabc,wabc,f0,q0,F_app,ax,ay,az)
global A E Ks G I sncon

K_sys=zeros(3*sNode_q,3*sNode_q);    F_sys=zeros(3*sNode_q,1);

for el=1:sncon
    n1=ncon(el,1);    n2=ncon(el,2);
    xa=Node(n1,1);    xc=Node(n2,1);
    
    wabc_e=wabc(el,:);     uabc_e=uabc(el,:);
    
    % K_ele=calculate_K_element(xa,xb,A,E,Ks,G,I,wab_e,ax,ay,az);
    [K_ele]=calculate_K__element(xa,xc,A,E,Ks,G,I,wabc_e,uabc_e,ax,ay,az);
    F_ele=calculate_F_element(xa,xc,f0,q0,F_app);
    
    % Making K_system and F_system
    r1=3*(n1-1+el-1)+1:3*(n1-1+el-1)+3;   % Note: Just for Unbranched frams
    r2=3*(n2-1+el-1)+1:3*(n2-1+el-1)+3;
    r3=3*(n2-1+el)+1:3*(n2-1+el)+3;
    
    rrr=[r1 r2 r3];
   
    K_sys(rrr(:),rrr(:))=K_sys(rrr(:),rrr(:))+K_ele(:,:);
    F_sys(rrr(:))=F_sys(rrr(:))+F_ele(:);
    
    
end

%% function  Assembling_KT_R
function [T_sys,R_sys]=T_R_sys(x_ele,sNode_q,Node,ncon,uabc,wabc,f0,q0,F_app,ax,ay,az)
global A E Ks G I sncon

T_sys=zeros(3*sNode_q,3*sNode_q);    F_sys=zeros(3*sNode_q,1);
R_sys=zeros(3*sNode_q,1);
for el=1:sncon
    n1=ncon(el,1);    n2=ncon(el,2);
    xa=Node(n1,1);    xc=Node(n2,1);
    
    wabc_e=wabc(el,:);     uabc_e=uabc(el,:);
    
    [K_ele T_ele]=calculate_K_T_element(xa,xc,A,E,Ks,G,I,wabc_e,uabc_e);
    F_ele=calculate_F_element(xa,xc,f0,q0,F_app);
    X_ele=x_ele(el,:);
    R_ele=-(F_ele-K_ele*X_ele');
    
    % Making K_system and F_system
    r1=3*(n1-1+el-1)+1:3*(n1-1+el-1)+3;   % Note: Just for Unbranched frams
    r2=3*(n2-1+el-1)+1:3*(n2-1+el-1)+3;
    r3=3*(n2-1+el)+1:3*(n2-1+el)+3;
    
    rrr=[r1 r2 r3];
    T_sys(rrr(:),rrr(:))=T_sys(rrr(:),rrr(:))+T_ele(:,:);
    R_sys(rrr(:))=R_sys(rrr(:))+R_ele(:);
    
end
h=1;

%% function  Assembling_KT_R
function [T_sys,R_sys]=T_R_sys_Dynamic(x_ele,sNode_q,Node,ncon,uabc,wabc,f0,q0,F_app,ax,ay,az,M_sys_con1,coef_damp,X,Xd,Xdd,XI,XdI,XddI,R_sys_con1)
global A E Ks G I sncon dt am ro_t
persistent gamma_ beta_
if isempty(beta_)
    gamma_=1/2;              beta_=1/4;
end
% [M_sys]=Mass_sys1(sNode_q,Node,ncon,am);
c0 = 1/(beta_*dt*dt) ;
c1 = gamma_/(beta_*dt) ;
c2 = 1/(beta_*dt) ;
c3 = 1/(beta_*2) - 1 ;
c4 = gamma_/beta_ - 1 ;
% c5 = 0.5*dt*(gamma_/beta_ - 2 ) ;
c5 = dt*(gamma_/beta_ - 1 ) ;
c6 = dt*(1 - gamma_ ) ;
c7 = dt* gamma_  ;
T_sys=zeros(3*sNode_q,3*sNode_q);    F_sys=zeros(3*sNode_q,1);
R_sys=zeros(3*sNode_q,1);
for el=1:sncon
    n1=ncon(el,1);    n2=ncon(el,2);
    xa=Node(n1,1);    xc=Node(n2,1);
    
    wabc_e=wabc(el,:);     uabc_e=uabc(el,:);
    
    % K_ele=calculate_K_element(xa,xb,A,E,Ks,G,I,wab_e,ax,ay,az);
    [K_ele T_ele]=calculate_K_T_element(xa,xc,A,E,Ks,G,I,wabc_e,uabc_e);
%     [K_ele]=calculate_K__element(xa,xc,A,E,Ks,G,I,wabc_e,uabc_e,ax,ay,az);
    F_ele=calculate_F_element(xa,xc,f0,q0,F_app);
    X_ele=x_ele(el,:);
%     R_ele=-(F_ele-K_ele*X_ele');
    
    % Making K_system and F_system
    r1=3*(n1-1+el-1)+1:3*(n1-1+el-1)+3;   % Note: Just for Unbranched frams
    r2=3*(n2-1+el-1)+1:3*(n2-1+el-1)+3;
    r3=3*(n2-1+el)+1:3*(n2-1+el)+3;

    rrr=[r1 r2 r3];
    
    Xii=XI(rrr);
    Xdii=XdI(rrr);
    Xddii=XddI(rrr);
    Xi=X(rrr);
    Xdi=Xd(rrr);
    Xddi=Xdd(rrr);
    R_1=R_sys_con1(rrr);
%     M=M_sys(rrr,rrr);
    [M_ele]=calculate_M_element(xa,xc,ro_t,A,I);
    M=M_ele*am;
    C=coef_damp*K_ele;
    Keff = c0*M + c1*C + K_ele ;
    Feff = F_ele+ M*(c0*Xi+c2*Xdi+c3*Xddi) ...
    +C*(c1*Xi+c4*Xdi+c5*Xddi) ;
%     R_ele=(Feff-Keff*X_ele');
%     R_ele=-(M*Xddi + C*Xdi + K_ele*X_ele')+F_ele;
%     R_ele=(M*Xddi + C*Xdi - Keff*X_ele')-Feff;
%     R_ele=-(F_ele-K_ele*X_ele');
%     R_ele=-(F_ele - Keff*X_ele' - M*Xddi - C*Xdi);
%     R_ele=-(F_ele - Keff*X_ele');
%     R_ele=(Feff - K_ele*X_ele');
%     R_ele=F_ele - M*Xddi - C*Xdi - K_ele*Xi;
%     R_ele=Feff - M*Xddi - C*Xdi - Keff*Xi;
%     DYNAMIC=M*(Xddi) + C*(Xdi);
%     DYNAMIC=M*(c3*Xddi) + C*(c4*Xdi);
    DYNAMIC=M*(c0*(Xii-Xi) - c2*Xdi - c3*Xddi) + C*(c1*(Xii-Xi) - c4*Xdi - c5*Xddi);
%     DYNAMIC=C*(c1*(Xii-Xi) - c4*Xdi - c5*Xddi);
%     P= M*Xddii + C*Xdii + Keff*Xii;
%     P= M*Xddi + C*Xdi + K_ele*(Xi);
    P= Keff*(Xii) + DYNAMIC;
%     P= K_ele*(Xi);
%     R_ele= F_ele - DYNAMIC - P;
    R_ele= -F_ele - P;
    l=4;
%     R_ele = F_ele - M*(c0*Xi+c2*Xdi+c3*Xddi) ...
%     - C*(c1*Xi+c4*Xdi+c5*Xddi) -K_ele*Xi;
%     R_ele= Feff - M*Xddi + C*Xdi + Keff*Xi;

%     R_ele=R_ele - R_1;
%     M*(c0*(Xi-Xii) - c2*Xdi - c3*Xddi)
%     C*(c1*(Xi-Xii) - c4*Xdi -c5*Xddi) 
%     K_ele*Xi
    T_sys(rrr(:),rrr(:))=T_sys(rrr(:),rrr(:))+T_ele(:,:);
    R_sys(rrr(:))=R_sys(rrr(:))+R_ele(:);
    
%     K_sys(rrr(:),rrr(:))=K_sys(rrr(:),rrr(:))+K_ele(:,:);
%     F_sys(rrr(:))=F_sys(rrr(:))+F_ele(:);
    
    
end
h=1;

%% function  Assembling_M
function [M_sys]=Mass_sys1(sNode_q,Node,ncon,am)
global A  I sncon ro_t

M_sys=zeros(3*sNode_q,3*sNode_q);

for el=1:sncon
    n1=ncon(el,1);    n2=ncon(el,2);
    xa=Node(n1,1);    xc=Node(n2,1);
    
    % Making K_system and F_system
    r1=3*(n1-1+el-1)+1:3*(n1-1+el-1)+3;   % Note: Just for Unbranched frams
    r2=3*(n2-1+el-1)+1:3*(n2-1+el-1)+3;
    r3=3*(n2-1+el)+1:3*(n2-1+el)+3;    
    rrr=[r1 r2 r3];
    
    [M_ele]=calculate_M_element(xa,xc,ro_t,A,I);    
    M_sys(rrr(:),rrr(:))=M_sys(rrr(:),rrr(:))+(am*M_ele(:,:));
end

%% function  Assembling_M
function [M_sys]=Mass_sys2(sNode_q,Node,ncon,am)
global A  I sncon ro_t dx

M_sys=zeros(3*sNode_q,3*sNode_q);
for el=1:sncon
    n1=ncon(el,1);    n2=ncon(el,2);
    xa=Node(n1,1);    xc=Node(n2,1);
    
    % Making K_system and F_system
    r1=3*(n1-1+el-1)+1:3*(n1-1+el-1)+3;   % Note: Just for Unbranched frams
    r2=3*(n2-1+el-1)+1:3*(n2-1+el-1)+3;
    r3=3*(n2-1+el)+1:3*(n2-1+el)+3;
    rrr=[r1 r2 r3];
    
    [M_ele]=mass_mms(xa,xc,ro_t,A,I,dx);
    M_sys(rrr(:),rrr(:))=M_sys(rrr(:),rrr(:))+(am*M_ele(:,:));
end

 %% function constraints Boundary K AND F
    function [K_sys_con1,F_sys_con1]=K_F_sys_con1(K_sys,F_sys,np_deg,sNode_q)  % apply constraints Boundary conds)
K_sys_con1=K_sys;
F_sys_con1=F_sys;
K_sys_con1(np_deg,:)=0;
K_sys_con1(:,np_deg)=0;
np_z_diag=zeros(3*sNode_q,1);
np_z_diag(np_deg)=1;
np_z_diag_1=diag(np_z_diag);
K_sys_con1=K_sys_con1+np_z_diag_1;
F_sys_con1(np_deg)=0;
 
%% function constraints Boundary KT AND R
    function [T_sys_con1,R_sys_con1]=KT_R_sys_con2(T_sys,R_sys,np_deg,sNode_q)  % apply constraints Boundary conds)
T_sys_con1=T_sys;
R_sys_con1=R_sys;
T_sys_con1(np_deg,:)=0;
T_sys_con1(:,np_deg)=0;
np_z_diag=zeros(3*sNode_q,1);
np_z_diag(np_deg)=1;
np_z_diag_1=diag(np_z_diag);
T_sys_con1=T_sys_con1+np_z_diag_1;
R_sys_con1(np_deg)=0;

%% function constraints Boundary Mass matrix
function [M_sys_con1]=Mass_sys_con1(M_sys,np_deg,sNode_q,am)  % apply constraints Boundary conds)
M_sys_con1=M_sys;
M_sys_con1(np_deg,:)=0;
M_sys_con1(:,np_deg)=0;
np_z_diag=zeros(3*sNode_q,1);
np_z_diag(np_deg)=1;
np_z_diag_1=diag(np_z_diag);
M_sys_con1=M_sys_con1+np_z_diag_1;
M_sys_con1=M_sys_con1*am;

 %% function static anlysis 
function [] =static(F_app,x_ele,f0,q0,X_load_incre,NI_F0,np_deg,load_increment_F0)
global  sncon ncon Node sNode_q ax ay az
wabc=zeros(sncon,3);      uabc=zeros(sncon,3);
% ipmact load Ray Clough
Load=zeros(1,NI_F0);
Load(1,1)=0;
Load(1,2)=5;
Load(1,3)=8;
Load(1,4)=7;
Load(1,5)=5;
Load(1,6)=3;
Load(1,7)=2;
Load(1,8)=1;
Load(1,9)=0;
for nom_inc_F0=1:NI_F0
%% The main while loop
% F_app=F_app+F0;
 F_app=Load(nom_inc_F0);
toler=100;   toler_w=100;   iter=0;
while toler>0.001 && iter<100
    iter=iter+1;
%% Assembling
[T_sys,R_sys]=T_R_sys(x_ele,sNode_q,Node,ncon,uabc,wabc,f0,q0,F_app,ax,ay,az);
%% Applying Boundary Conditions
[T_sys_con1,R_sys_con1]=KT_R_sys_con2(T_sys,R_sys,np_deg,sNode_q);  % apply constraints Boundary conds)
Determinant(nom_inc_F0)=det(1e-6*T_sys_con1);
% ----------------

dx_sys=-T_sys_con1\R_sys_con1;

if iter==1 && nom_inc_F0==1
X_sys=dx_sys;
elseif iter==1
X_sys=X_load_incre(:,nom_inc_F0-1)+dx_sys;
else
X_sys=X_sys_former+dx_sys;
end    %if

toler=sum(abs(R_sys_con1));    %          toler=max(abs(dx_sys));
X_sys_former=X_sys;

%% calculation of wa, wb which are the only nonlinear parameters
for el=1:sncon
  n1=ncon(el,1);
  n2=ncon(el,2);
%   r1=3*(n1-1)+1:3*(n1-1)+3;
%   r2=3*(n2-1)+1:3*(n2-1)+3;
%   rr=[r1 r2];   

  r1=3*(n1-1+el-1)+1:3*(n1-1+el-1)+3;   % Note: Just for Unbranched frams
  r2=3*(n2-1+el-1)+1:3*(n2-1+el-1)+3;
  r3=3*(n2-1+el)+1:3*(n2-1+el)+3;
 
  rr=[r1 r2 r3];  
  x_ele(el,:)=X_sys(rr);
  wabc(el,:)=x_ele(el,[2 5 8]);
  uabc(el,:)=x_ele(el,[1 4 7]);
end

end     %while

    X_load_incre(:,nom_inc_F0)=X_sys;

% inp_graph(nom_inc_F0)=X_load_incre((sNode_q-1)/2*3+1,nom_inc_F0);
inp_graph(nom_inc_F0)=X_load_incre(1,nom_inc_F0);
end     %for nom_inc_F0
%% Showing the element displacements for the maximum applied load.
NI=NI_F0;
hold on
plot(inp_graph,load_increment_F0*(0:NI_F0-1))

 %% function Newmark 
function [X,Xd,Xdd] =Newmark(M,C,K,F,dt,Xi,Xdi,Xddi,sNode_q,n)
persistent gamma_ beta_
if isempty(beta_)
    gamma_=1/2;              beta_=1/4;
end
c0 = 1/(beta_*dt*dt) ;
c1 = gamma_/(beta_*dt) ;
c2 = 1/(beta_*dt) ;
c3 = 1/(beta_*2) - 1 ;
c4 = gamma_/beta_ - 1 ;
c5 = 0.5*dt*(gamma_/beta_ - 2 ) ;
% c5 = dt*(gamma_/beta_ - 1 ) ;
c6 = dt*(1 - gamma_ ) ;
c7 = dt* gamma_  ;
Keff = c0*M + c1*C + K ;
% Feff = F+ M*(c0*Xi+c2*Xdi+c3*Xddi) ...
%     +C*(c1*Xi+c4*Xdi+c5*Xddi) ;
Feff = F;

%  Mii=M((sNode_q-1)/2,(sNode_q-1)/2);
% if Mii==0
%     X=Keff\Feff;
%     X=X+Xi;
% else
    X=Keff\Feff;
    X=X+Xi;
% end

Xdd= c0*(X-Xi) - c2*Xdi - c3*Xddi ;
Xd = Xdi + c6*Xddi + c7*Xdd ;



%% function dynamic 
function [] =dynamic(nt,v0,coef_damp,np_deg,max_iter,F0,x_ele,f0,q0)
    global  sncon ncon Node sNode_q ax ay az am dt Ty
    % ipmact load Ray Clough
    Load=zeros(1,nt);
    Load(1,1)=0;
    Load(1,2)=5;
    Load(1,3)=8;
    Load(1,4)=7;
    Load(1,5)=5;
    Load(1,6)=3;
    Load(1,7)=2;
    Load(1,8)=1;
    Load(1,9)=0;
    Xmat=zeros(3*sNode_q,nt);
    Xdmat=zeros(3*sNode_q,nt);
    Xddmat=zeros(3*sNode_q,nt);
%     Xmat(1,1)=0.2;
%     Xdmat(1,1)=5;
% Assembling mass matrix
    [M_sys]=Mass_sys1(sNode_q,Node,ncon,am);
%     [M_sys]=Mass_sys2(sNode_q,Node,ncon,am);
% Applying Boundary Conditions mass matrix
    [M_sys_con1]=Mass_sys_con1(M_sys,np_deg,sNode_q,am);
    for n=1:nt-1
%          F_app=10*sin(5*n*dt/Ty);
%          F_app=10;
%          F_app=0;
         F_app=Load(n+1);
%          F_app=F_app+F0;
        toler=100;   toler_w=100;   iter=0;
        if n==1
            X_sys=zeros(3*sNode_q,1);
%             [X_sys]=Boundary(E,ro_t,L,v0,sncon,dx,sNode_q,r,stress_cr);
            Xmat(:,1)=X_sys;
            dXmat(1,1)=v0;
            [wabc,uabc,x_ele]=uw_ele(ncon,X_sys,sncon);
            R_sys_con1=zeros(3*sNode_q,1);
        end %end if n==1

        while  toler>0.001 && iter<max_iter
            iter=iter+1;
            
%         Assembling
            Xi=Xmat(:,n);
            Xdi=Xdmat(:,n);
            Xddi=Xddmat(:,n);
            XI=Xmat(:,n+1);
            XdI=Xdmat(:,n+1);
            XddI=Xddmat(:,n+1);
            
%          Mii=M_sys_con1((sNode_q-1)/2,(sNode_q-1)/2);%dynamic analysis
%             if Mii==0% static analysis
%                 XI=Xmat(:,n+1);
%                 XdI=Xdmat(:,n+1);
%                 XddI=Xddmat(:,n+1);
%             else% dynamic analysis
                if iter==1
                    XI=Xmat(:,n);
                    XdI=Xdmat(:,n);
                    XddI=Xddmat(:,n);
%                 else
%                     XI=Xmat(:,n+1);
%                     XdI=Xdmat(:,n+1);
%                     XddI=Xddmat(:,n+1);
                end
%             end
            
            if n==1 
                XI=zeros(3*sNode_q,1);
            else
            XI=Xmat(:,n-1);
            end
            
            [T_sys,R_sys]=T_R_sys_Dynamic(x_ele,sNode_q,Node,ncon,uabc,wabc,f0,q0,F_app,ax,ay,az,M_sys_con1,coef_damp,Xi,Xdi,Xddi,XI,XdI,XddI,R_sys_con1);
%             [T_sys,R_sys]=T_R_sys(x_ele,sNode_q,Node,ncon,uabc,wabc,f0,q0,F_app,ax,ay,az)%         Applying Boundary Conditions
            [T_sys_con1,R_sys_con1]=KT_R_sys_con2(T_sys,R_sys,np_deg,sNode_q);  % apply constraints Boundary conds)
            C=coef_damp*T_sys_con1;
             X_sys_former=X_sys;     
%              f=F_sys_con1;
%              k3=K_sys_con1;
             f=R_sys_con1;
             k3=T_sys_con1;
            
            if n==1 && iter==1;
                  dddXmat=(f-k3*Xi-C*Xdi)\M_sys_con1;
                  ddXmat(:,1)=dddXmat';
            end
            [X,Xd,Xdd] =Newmark(M_sys_con1,C,k3,f,dt,Xi,Xdi,Xddi,sNode_q,n);
            
            XI=X;
            XdI=Xd;
            XddI=Xdd;
            
            X_sys=X;
%             [wabc,uabc,x_ele]=uw_ele(ncon,X_sys,sncon);
            for el=1:sncon
          n1=ncon(el,1);
          n2=ncon(el,2);
%           r1=3*(n1-1)+1:3*(n1-1)+3;
%           r2=3*(n2-1)+1:3*(n2-1)+3;
%           rr=[r1 r2];  

          r1=3*(n1-1+el-1)+1:3*(n1-1+el-1)+3;   % Note: Just for Unbranched frams
          r2=3*(n2-1+el-1)+1:3*(n2-1+el-1)+3;
          r3=3*(n2-1+el)+1:3*(n2-1+el)+3;

          rr=[r1 r2 r3];  
          x_ele(el,:)=X_sys(rr);

          wabc(el,:)=x_ele(el,[2 5 8]);

          uabc(el,:)=x_ele(el,[1 4 7]);
        end
%             contorol
            toler=max(abs((X_sys_former-X_sys)./X_sys));

            if iter>max_iter
                break
            end

        end
        Xmat(:,n+1)=X;
        Xdmat(:,n+1)=Xd;
        Xddmat(:,n+1)=Xdd;
    end
        
hold on
draw= Xmat((sNode_q-1)/2*3+1,:);
draw2= Xmat((sNode_q-1)/2*3+2,:);
draw3= Xmat((sNode_q-1)/2*3+3,:);
draw1= Xmat(1,:);
draw12= Xmat(2,:);
draw13= Xmat(3,:);
% draw= Xmat(1,:);
% draw1= E*(Xmat(sNode_q*3-2,:)-Xmat(sNode_q*3-8,:))/dx;
% draw2= E*(Xmat(3*3-2,:)-Xmat(3*3-8,:))/dx;
% tvec=0:dt:TT;
% plot(tvec,draw);
plot(draw1);
% plot(draw3);

% plot(draw2);
% plot(draw12);
% plot(draw3);
% plot(draw13);


%% Gauss integration 
 function Gauss=gauss(NGP)

% gauss=cell(5,1);
% gauss{1}=[0 2];
% gauss{2}=[0.5773502692 1
%                -0.5773502692 1];
% gauss{3}=[0.7745966692 0.5555555555
%                 0 0.8888888889
%                 -0.7745966692 0.5555555555];
% gauss{4}=[0.3399810435 0.6521451548
%                 0.8611363116 0.3478548451
%                 -0.3399810435 0.6521451548
%                 -0.8611363116 0.3478548451];
% gauss{5}=[0.5384693101 0.4786286705
%                 0.9061798459 0.2369268850
%                 0.0000000000 0.5688888889
%                 -0.5384693101 0.4786286705
%                 -0.9061798459 0.2369268850];    
            
gauss{1}=[0 2];
gauss{2}=[1/sqrt(3) 1
               -1/sqrt(3) 1];
gauss{3}=[sqrt(0.6) 5/9
                0 8/9
                -sqrt(0.6) 5/9];
            const1=sqrt(4.8);
gauss{4}=[sqrt((3+const1)/7) 0.5-1/(3*const1)
                sqrt((3-const1)/7) 0.5+1/(3*const1)
                -sqrt((3-const1)/7) 0.5+1/(3*const1)
                -sqrt((3+const1)/7) 0.5-1/(3*const1)];
gauss{5}=[0.5384693101 0.4786286705
                0.9061798459 0.2369268850
                0.0000000000 0.5688888889
                -0.5384693101 0.4786286705
                -0.9061798459 0.2369268850];     
   Gauss=gauss{NGP};                   
 return

%% function input6
     function [Node ncon Constraints]=input6()
c1=1e270;
c2=1e-270;
%[node_number rigidity_u rigidity_w rigidity_teta] 

Constraints=[1	c2	c1	c2
21	c1	c1	c2];

force=[];

%%Node  &   ncon
poisson_ratio=0 ; %No Pois. Coupl. poisson_ratio=0    Poi. Coup. included    =0.3

%%ncon
%3: inner diameter 4:epsilon=raphness   5:E  6:G   7:Rho   8:thickness  9:poisson ratio
Node=[								
0	0	0
5	0	0
10	0	0
15	0	0
20	0	0
25	0	0
30	0	0
35	0	0
40	0	0
45	0	0
50	0	0
55	0	0
60	0	0
65	0	0
70	0	0
75	0	0
80	0	0
85	0	0
90	0	0
95	0	0
100	0	0
];								
ncon=[								
1	2	0.797	0.00004	30000000	12000000	7900	0.008	0.3
2	3	0.797	0.00004	30000000	12000000	7900	0.008	0.3
3	4	0.797	0.00004	30000000	12000000	7900	0.008	0.3
4	5	0.797	0.00004	30000000	12000000	7900	0.008	0.3
5	6	0.797	0.00004	30000000	12000000	7900	0.008	0.3
6	7	0.797	0.00004	30000000	12000000	7900	0.008	0.3
7	8	0.797	0.00004	30000000	12000000	7900	0.008	0.3
8	9	0.797	0.00004	30000000	12000000	7900	0.008	0.3
9	10	0.797	0.00004	30000000	12000000	7900	0.008	0.3
10	11	0.797	0.00004	30000000	12000000	7900	0.008	0.3
11	12	0.797	0.00004	30000000	12000000	7900	0.008	0.3
12	13	0.797	0.00004	30000000	12000000	7900	0.008	0.3
13	14	0.797	0.00004	30000000	12000000	7900	0.008	0.3
14	15	0.797	0.00004	30000000	12000000	7900	0.008	0.3
15	16	0.797	0.00004	30000000	12000000	7900	0.008	0.3
16	17	0.797	0.00004	30000000	12000000	7900	0.008	0.3
17	18	0.797	0.00004	30000000	12000000	7900	0.008	0.3
18	19	0.797	0.00004	30000000	12000000	7900	0.008	0.3
19	20	0.797	0.00004	30000000	12000000	7900	0.008	0.3
20	21	0.797	0.00004	30000000	12000000	7900	0.008	0.3
];								

 return 

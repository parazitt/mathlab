%   IN THE NAME OF GOD
%   The method Newton Rophsun and famley Newmark for displesmen
%   The method finite elements for copclitment
%   The method uoler bernoly beam for matrix mass and stiffness
%   wighit mohamad mahdi salati AND github
%   fo perblem impact and dynamic 

% investigation for axial loads
%    The Plane Frame Element K M  Global Matrises
%in the pro for stress impact
%stress for impact mass to bar and impact mass to bar free and impact bar
%to bar.
% newmark metod and matrex mass and stiffness of uoler bernoly beam
% matrex mass and stiffness is coefficients of parameter
% dar in try mami khahim ba amale sefr nemodan matrix sakhti ke moj tanesh
% be an vared nashodeh ast  nataij ra be dast avarim  v faqat sorat moase
% khahad bod.


function []=static_nonlinear()
clc
[Node ncon Constraints]= input6();
%% Preliminaries of the main while loop

global A E Ks G I sncon ro_t nsec dx sNode_q
Ks=5/6;
Ks=1; % Ks is applied in ay and az.

%% deta 1
E= 29e9  ;
E= 5*720000/11  ;
% ro_t= 2500;
ro_t= 1/9;
nu=0.2;
bbbb=0.5;
hhhh=0.3;
L=6;
M2=ro_t*bbbb*hhhh*L;
M1=45*M2;

%% deta 2
% E= 210e9  ;
% ro_t=7800;
% stress_cr=1500e6;
% L=0.193;
% nu=0.2;
% bbbb=8.7*10^-3;
% hhhh=0.63*10^-3;
% M2=ro_t*bbbb*hhhh*L;
% M1=0.369;
%% calc
A=bbbb*hhhh;
ay=1.22;
az=1.2;
Ay=(1/ay)*A;
Az=(1/az)*A;
A=Ay;
Iz=bbbb*hhhh^3/12;
K_kol=3*E*Iz/L^3;
ww=sqrt(K_kol/M2);
I=Iz;
r=sqrt(I/A);
ct=sqrt(E/ro_t);
nx=nsec+1;
dx=L/nsec;
% dt=dx/ct;
dt=0.1;
% Ty=L/ct;
Ty=1;
TT=1*Ty;
% dt=dx/10/ct;
% coef_damp=0.000003;
coef_damp=0;
% coef_damp=0.01;
nt=floor(TT/dt)+1;
nty=floor(Ty/dt)+1;

Pcr=10*pi*pi*E*Iz/(L*L)     ;
G=E/(2*(1+nu));

max_iter=50;
%% matrix stiffness Euler Beam
Axx=E*A;
Dxx=E*I;
Sxx=Ks*G*A;
%coefficients of parameter (ax,ay,az,=>stiffness matrix and am=>mass matrix)
ax=1;   ay=1;   az=1;   am=0;

%% 
% v0=4;
% v0=6.5;
% v0=15;
v0=0;
% v0=17;
% v0=19;%cretical
% v0=18;
Load=zeros(1,nt);
Load(1,1)=0
Load(1,2)=5;
Load(1,3)=8;
Load(1,4)=7;
Load(1,5)=5;
Load(1,6)=3;
Load(1,7)=2;
Load(1,8)=1;
Load(1,9)=0;
% Load(:,:)=0.001;
%% Axial impact of a rigid mass to the free end of a long rod cantilever
% [Load]=mass_bar_free(L,ct,nt,dt,Ty,A,v0,E,ro_t,M1,M2);
% [v]=mass_bar_freev(L,ct,nt,dt,Ty,v0,M1,M2);
%% Axial impact of a rigid mass to a free end of the long rod free ends
% [Load]=mass_bar(L,ct,nt,dt,Ty,A,v0,E,ro_t,M1,M2);
% [v]=mass_barv(L,ct,nt,dt,Ty,v0,M1,M2);
%% Axial impact of a subjected to a free end of the long rod cantilever
% [Load]=bar_bar(L,ct,nt,dt,Ty,A,v0,E,ro_t,M1,M2);
% [v]=bar_barv(L,ct,nt,dt,Ty,v0,M1,M2);
% Load(:,:)=100;%test dynamic
%% deta matrix
sncon=size(ncon,1);
sNode=size(Node,1);   sNode_q=2*sNode-1;
wabc=zeros(sncon,3);      uabc=zeros(sncon,3);
% K_sys=zeros(3*sNode,3*sNode);    F_sys=zeros(3*sNode,1);
Node=zeros(sNode,3);
harp=0:dx:L;    harp=harp';
Node(:,1)=harp;


X_sys_former=zeros(3*sNode_q,1);
x_ele=zeros(sncon,9);  %10*ones(sncon,6);
%% Applying Boundary Conditions
c1=1e27;
c2=1e-27;
Constraints=[1	c2	c2	c2
    nsec+1	c1	c1	c1];

[as bs]=find(Constraints>1000);
np_deg=sort(2*3*(-1+Constraints(as))+bs-ones(size(bs,1),1));

%% LOADING
f0=0;
% q=0.0001;
q=0;
load_increment_q0=0.0001;
% load_increment_q0=0;
BUCKLING=2466;
% F0_=3500;
F0_=BUCKLING;
stress_cr=Pcr/A;
load_increment_F0=1;
% F0=0.01;
% load_increment_F0=0.01;

NI_F0=F0_/load_increment_F0;
inp_graph=zeros(NI_F0,1);
F0=load_increment_F0;
Determinant=zeros(NI_F0,1);

NI_q0=q/load_increment_q0;
q0=load_increment_q0;
X_load_incre=zeros(3*sNode_q,NI_F0);
F_app=0;

Xmat=zeros(3*sNode_q,nt);
Xdmat=zeros(3*sNode_q,nt);
Xddmat=zeros(3*sNode_q,nt);
Xdmat(1,1)=v0;
check=zeros(nt,2);
%% mass matrix
[M_sys]=Mass_sys1(sNode_q,Node,ncon,am);
% [M_sys]=Mass_sys2(sNode_q,Node,ncon,am);
% [M_sys]=Mass_sys3(sNode_q,Node,ncon,am,dx);
[M_sys_con1]=Mass_sys_con1(M_sys,np_deg,sNode_q) ; % apply constraints Boundary conds)
[M_sys_con1]=am*[M_sys_con1];



results=zeros(nt,1);
ii=1;
jj=1;
for n=1:nt
    %% The main while loop
    %         F_app=F_app+F0;   % F_app: applied load at the current load incerement
    F_app=F_app+Load(n);%k dynamic test step load
%     F_app=Load(n+1);
%     F_app=0.001;
%         detlaF=F_app-Load(n);

%% IMPACT LOAD
%     delta=n*dt-dt;
% %     hex=101*dt;
% hex=0.025;
%     if  delta<=hex
%         F_app=3170300*(delta/hex);
%     elseif delta<=2*hex
%         F_app=3170300*(1-((delta-hex)/hex));
%     else
%         F_app=0;
%     end


    toler=100;   toler_w=100;   iter=0;
    X_sys=zeros(3*sNode_q,1);

    if n==1

        dXmat(1,1)=v0;
        X_sys=zeros(3*sNode_q,1);
%         [X_sys]=Boundary(E,ro_t,L,v0,sncon,dx,sNode_q,r,stress_cr);
%         Xmat(:,n)=X_sys;
%         Xmat(1,1)=v0*dt;

        [wabc,uabc]=uw_ele(ncon,X_sys);
        [F]=F_M_V1(X_sys,ncon);

    else
        [F]=F_M_V1(X,ncon);
    end
%     if rem(ii,2*sncon)==0
%         ii=1;
%         jj=1;
%     end
%     if ii<=sncon
%         jj=jj+6;
%     else
%         jj=jj-6;
%     end
%     if n==1
%         dXmat(1,1)=v0;
%     else
%         dXmat(jj,n)=v(n+1);
%     end

    while toler>0.001 && iter<max_iter
        iter=iter+1;
        %% Assembling
        [K_sys,F_sys]=K_F_sys(sNode_q,Node,ncon,uabc,wabc,f0,q0,F_app,ax,ay,az);%non load +
%         [K_sys,F_sys]=K_F_sys1(sNode_q,Node,ncon,uabc,wabc,f0,q0,F_app,ax,ay,az);%LOAD +
%         [K_sys,F_sys]=K_F_sys2(sNode_q,Node,ncon,uabc,wabc,f0,q0,F_app,ax,ay,az,F);%F_M_V+
        
        [K_sys_con1,F_sys_con1]=K_F_sys_con1(K_sys,F_sys,np_deg,sNode_q);  % apply constraints Boundary conds)
        C=coef_damp*K_sys_con1;
        X_sys_former=X_sys;
        
        f=F_sys_con1;
        k3=K_sys_con1;

        Xi=Xmat(:,n);
        Xdi=Xdmat(:,n);
        Xddi=Xddmat(:,n);
        [X,Xd,Xdd] =Newmark(M_sys_con1,C,k3,f,dt,Xi,Xdi,Xddi,sNode_q,n,nty);
        Xmat(:,n+1)=X;
        Xdmat(:,n+1)=Xd;
        Xddmat(:,n+1)=Xdd;
        
        X_sys=X;
        toler=max(abs((X_sys_former-X_sys)./X_sys));
        
        %% calculation of wa, wb which are the only nonlinear parameters
        [wabc,uabc]=uw_ele(ncon,X_sys);
       
        if iter==max_iter
            break
        end
    end     %while
    
    Xmat(:,n+1)=X_sys;
    results(n+1)=iter;

end
%% plot
% inp_graph(nom_inc_F0)=X_load_incre((sNode_q-1)/2*3+2,n+1);
% inp_graph(nom_inc_F0)=X_load_incre((sNode_q-1)/2*3+2,nom_inc_F0);


% if nom_inc_F0==1
%     X_load_incre(:,nom_inc_F0)=X_sys;
% else
% X_load_incre(:,nom_inc_F0)=X_load_incre(:,nom_inc_F0-1)+X_sys;
% end    %if
%     Xmat(:,nt)=X_sys;
%     [wabc,uabc]=uw_ele(ncon,X_sys);
%       hold on
% plot(uabc);

     
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
plot(draw);
plot(draw1);
plot(draw2);
plot(draw12);
plot(draw3);
plot(draw13);

% plot(tvec,draw1);
% plot(tvec,draw2);
% plot(tvec,Load);
return

function [F]=F_M_V1(X_sys,ncon);
global sncon E I A dx sNode_q
% wabc
% uabc
    F=zeros(3*sNode_q,1);
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
    du=(x_ele(el,1)-x_ele(el,7))/dx;
    dw=(x_ele(el,2)-x_ele(el,8))/dx;
    dww=(x_ele(el,2)-x_ele(el,8))/dx^2;
    N=E*A*(du+(0.5*dw^2));
    M=E*I*(dww/(1+dw^2));
    F_ele(1,1)=N;
    F_ele(2,1)=M;
    F_ele(9,1)=0;
    
    F(rr(:))=F(rr(:))+F_ele(:);
    %   xd_ele(el,:)=xd(rr);
    %   xdd_ele(el,:)=xdd(rr);
    %
    %   x_ele_loc(el,:)=RT(:,:,el)*x_ele(el,:)';
    %   xd_ele_loc(el,:)=RT(:,:,el)*xd_ele(el,:)';
    %   xdd_ele_loc(el,:)=RT(:,:,el)*xdd_ele(el,:)';
    wabc(el,:)=x_ele(el,[2 5 8]);
    
    
    %   if wab(el,2)<0.51771 && wab(el,2)>0.51769;
    %       wab(el,2)
    %   end
    uabc(el,:)=x_ele(el,[1 4 7]);
end
 l=4;   
return


function [Xmat]=Boundary(E,ro,L,v0,sncon,dx,sNode_q,r,stress_cr)

        %problem        
        
        ct=sqrt(E/ro);
        n=1;%replay stress wave
        tt=n*(L/ct);
        t=L/ct;

        mu=0.005;
        k1=10*3.1416^2;%loose clamp conditions
        % k1=5*3.1416^2;%simply supported conditions
        k2=3*3.1416^2;%loose clamp conditions
        % k2=2*3.1416^2;%simply supported conditions

        dx=L/sncon;
        dt=dx/ct;
        nt=tt/dt;
        lcr=r*sqrt(k1*E/stress_cr);
    %     lcr=36.66*10^-3;%loose clamp conditions
    %     lcr=14.46*10^-3;%simply supported conditions
        tcr=lcr/ct;
            for j=1:(sNode_q);
        
                r1=3*(j-1)+1;    r2=3*(j-1)+2;      r3=3*(j-1)+3;
                mmm=j-1;
                ttt= j*dt;
                xxx=mmm*dx;
                zeta=xxx/(lcr);
                TTT=ttt/tcr;
                www=r*k2/lcr;

                y=3*sin(pi*zeta)-sin(pi*3*zeta);%loose clamp conditions
%                 y=2*sin(pi*zeta)+sin(pi*2*zeta);%simply supported conditions
%                 Xmat(1,1)=v0*dt;

                    if xxx<=lcr
                        Xmat(r2,1)=mu*exp(www*(TTT-1))*y;
                    else
                        Xmat(r2,1)=0;
                    end
                        Xmat(r3,1)=0;

            end

function [Load]=mass_bar_free(L,ct,nt,dt,Ty,A,v0,E,ro_t,M1,M2)
for i=1:nt
    tt=i*dt;
    if tt<(Ty*2)
        AAA=-A*v0*sqrt(E*ro_t)*exp(-M2/M1*tt*ct/L);
    elseif tt<(Ty*4)
        ttt=tt-(Ty*2);
        AAA=-A*v0*sqrt(E*ro_t)*exp(-M2/M1*ttt*ct/L)*(2+exp(-2*M2/M1)-(2*ttt*M2*ct/(M1*L)));
%     elseif tt<(Ty*6)
%         if (M2/M1)>0.58;
%             if tt <(Ty*5)
%                 AAA=-A*v0*sqrt(E*ro_t)*exp(-M2/M1*2)*(2+exp(-2*M2/M1)-(2*2*M2/(M1)));
%             else
%                 AAA=2*(-A*v0*sqrt(E*ro_t)*exp(-M2/M1*2)*(2+exp(-2*M2/M1)-(2*2*M2/(M1))));
%             end
%         else
%             AAA=
    end
    Load(i)=AAA;
end
   

function [Load]=mass_bar(L,ct,nt,dt,Ty,A,v0,E,ro_t,M1,M2)
for i=1:nt
    tt=i*dt;
    if tt<=Ty
        AAA=-A*v0*sqrt(E*ro_t)*exp(-(M2/M1)*tt*(ct/L));
%         AAA_=-A*v0*sqrt(E*ro_t)*exp(-(M2/M1)*tt*(ct/L));
    elseif tt<=(Ty*2)
        ttt=tt-Ty;
        x_=ct*ttt;
%         x_bar=L-(x_/2);
        X_=L-x_;
        x_bar=X_;
        
%         AAA=A*v0*sqrt(E*ro_t)*exp(-(M2/M1)*tt*(ct/L));
%         AAA_=A*v0*sqrt(E*ro_t)*(exp(-M2*(x_bar-X_)/(M1*L))-exp(-M2/M1*(2-((x_bar+X_)/L))));%Stress at some point in the distance x from impactor
        AAA=A*v0*sqrt(E*ro_t)*(-1+exp(-M2/M1*(2-((x_bar+X_)/L))));%Stress at some point in the distance x from impactor
%         AAA=A*v0*sqrt(E*ro_t)*exp(-M2/M1*tt*ct/L);
%         AAA=-A*v0*sqrt(E*ro_t)*exp(-M2/M1*ttt*ct/L)*(2+exp(-2*M2/M1)-(2*ttt*M2*ct/(M1*L)));
%         AAA=A*v0*sqrt(E*ro_t)*exp(-(M2/M1)*(2-(X_/L)));
    elseif tt<=(3*Ty)
%         AAA=-A*v0*sqrt(E*ro_t)*(1-exp(-2*M2/M1));
        AAA=-A*v0*sqrt(E*ro_t)*(1-exp(-2*M2/M1));
%         AAA=-A*v0*sqrt(E*ro_t)*exp(-M2/M1*ttt*ct/L)*(exp(-2*M2/M1)-(2*ttt*M2*ct/(M1*L)));
    elseif tt<=(4*Ty)
        AAA=0;
    else
%         glk=A*v0*sqrt(E*ro_t)*exp(-2*M2/M1);
        AAA=0;
    end
    Load(i)=AAA;
end
% plot(Load)
% plot(Load_)
% glk=A*v0*sqrt(E*ro_t)*exp(-2*M2/M1);
% a=0

function [Load]=bar_bar(L,ct,nt,dt,Ty,A,v0,E,ro_t,M1,M2)
Q=0.5; 
% Q=1;
L2=Q*L;
for i=1:nt
    tt=i*dt;
    if tt<Ty
        AAA=-A*v0*sqrt(E*ro_t);
    elseif tt<=(Ty*2)
        ttt=tt-Ty;
        x_=ct*ttt;
%         x_bar=L-(x_/2);
        X_=L-x_;
        x_bar=X_;
        if Q==1
        AAA=-2*A*v0*sqrt(E*ro_t);
        else
            if (Q*L)>x_
                AAA=-2*A*v0*sqrt(E*ro_t);
            else
                AAA=-A*v0*sqrt(E*ro_t);
            end
        end
    elseif tt<=(Ty*3)
        ttt=tt-(Ty*2);
        x_=ct*ttt;
        if Q==1
        AAA=A*v0*sqrt(E*ro_t);
        else
            if (Q*L)>x_
                AAA=A*v0*sqrt(E*ro_t);
            else
                AAA=-A*v0*sqrt(E*ro_t);
            end
        end
    elseif tt<=(Ty*4)
        ttt=tt-(Ty*3);
        x_=ct*ttt;
        if Q==1
        AAA=2*A*v0*sqrt(E*ro_t);
        else
            if (Q*L)>x_
                AAA=-2*A*v0*sqrt(E*ro_t);
            else
                AAA=-A*v0*sqrt(E*ro_t);
            end
        end
    else
        AAA=0;
    end
    Load(i)=AAA;
end


function [v]=mass_bar_freev(L,ct,nt,dt,Ty,v0,M1,M2)
for i=1:nt
    tt=i*dt;
    if tt<(Ty*2)
        AAA=-v0*exp(-M2/M1*tt*ct/L);
    elseif tt<(Ty*4)
        ttt=tt-(Ty*2);
        AAA=-v0*exp(-M2/M1*ttt*ct/L)*(2+exp(-2*M2/M1)-(2*ttt*M2*ct/(M1*L)));
%     elseif tt<(Ty*6)
%         if (M2/M1)>0.58;
%             if tt <(Ty*5)
%                 AAA=-A*v0*sqrt(E*ro_t)*exp(-M2/M1*2)*(2+exp(-2*M2/M1)-(2*2*M2/(M1)));
%             else
%                 AAA=2*(-A*v0*sqrt(E*ro_t)*exp(-M2/M1*2)*(2+exp(-2*M2/M1)-(2*2*M2/(M1))));
%             end
%         else
%             AAA=
    end
    v(i)=AAA;
end
   

function [v]=mass_barv(L,ct,nt,dt,Ty,v0,M1,M2)
for i=1:nt
    tt=i*dt;
    if tt<=Ty
        AAA=-v0*exp(-(M2/M1)*tt*(ct/L));
%         AAA_=-A*v0*sqrt(E*ro_t)*exp(-(M2/M1)*tt*(ct/L));
    elseif tt<=(Ty*2)
        ttt=tt-Ty;
        x_=ct*ttt;
%         x_bar=L-(x_/2);
        X_=L-x_;
        x_bar=X_;
        
%         AAA=A*v0*sqrt(E*ro_t)*exp(-(M2/M1)*tt*(ct/L));
%         AAA_=A*v0*sqrt(E*ro_t)*(exp(-M2*(x_bar-X_)/(M1*L))-exp(-M2/M1*(2-((x_bar+X_)/L))));%Stress at some point in the distance x from impactor
        AAA=v0*(-1+exp(-M2/M1*(2-((x_bar+X_)/L))));%Stress at some point in the distance x from impactor
%         AAA=A*v0*sqrt(E*ro_t)*exp(-M2/M1*tt*ct/L);
%         AAA=-A*v0*sqrt(E*ro_t)*exp(-M2/M1*ttt*ct/L)*(2+exp(-2*M2/M1)-(2*ttt*M2*ct/(M1*L)));
%         AAA=A*v0*sqrt(E*ro_t)*exp(-(M2/M1)*(2-(X_/L)));
    elseif tt<=(3*Ty)
%         AAA=-A*v0*sqrt(E*ro_t)*(1-exp(-2*M2/M1));
        AAA=-v0*(1-exp(-2*M2/M1));
%         AAA=-A*v0*sqrt(E*ro_t)*exp(-M2/M1*ttt*ct/L)*(exp(-2*M2/M1)-(2*ttt*M2*ct/(M1*L)));
    elseif tt<=(4*Ty)
        AAA=0;
    else
%         glk=A*v0*sqrt(E*ro_t)*exp(-2*M2/M1);
        AAA=0;
    end
    v(i)=AAA;
end

function [v]=bar_barv(L,ct,nt,dt,Ty,v0,M1,M2)
Q=0.5; 
% Q=1;
L2=Q*L;
for i=1:nt
    tt=i*dt;
    if tt<Ty
        AAA=-v0;
    elseif tt<=(Ty*2)
        ttt=tt-Ty;
        x_=ct*ttt;
%         x_bar=L-(x_/2);
        X_=L-x_;
        x_bar=X_;
        if Q==1
        AAA=-2*v0;
        else
            if (Q*L)>x_
                AAA=-2*v0;
            else
                AAA=-v0;
            end
        end
    elseif tt<=(Ty*3)
        ttt=tt-(Ty*2);
        x_=ct*ttt;
        if Q==1
        AAA=v0;
        else
            if (Q*L)>x_
                AAA=v0;
            else
                AAA=-v0;
            end
        end
    elseif tt<=(Ty*4)
        ttt=tt-(Ty*3);
        x_=ct*ttt;
        if Q==1
        AAA=2*v0;
        else
            if (Q*L)>x_
                AAA=-2*v0;
            else
                AAA=-v0;
            end
        end
    else
        AAA=0;
    end
    v(i)=AAA;
end


function [wabc,uabc]=uw_ele(ncon,X_sys)  % apply constraints Boundary conds)
global sncon
% wabc
% uabc
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
    %   xd_ele(el,:)=xd(rr);
    %   xdd_ele(el,:)=xdd(rr);
    %
    %   x_ele_loc(el,:)=RT(:,:,el)*x_ele(el,:)';
    %   xd_ele_loc(el,:)=RT(:,:,el)*xd_ele(el,:)';
    %   xdd_ele_loc(el,:)=RT(:,:,el)*xdd_ele(el,:)';
    wabc(el,:)=x_ele(el,[2 5 8]);
    
    
    %   if wab(el,2)<0.51771 && wab(el,2)>0.51769;
    %       wab(el,2)
    %   end
    uabc(el,:)=x_ele(el,[1 4 7]);
end

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

function [M_sys_con1]=Mass_sys_con1(M_sys,np_deg,sNode_q)  % apply constraints Boundary conds)

M_sys_con1=M_sys;
% F_sys_con1=F_sys;
M_sys_con1(np_deg,:)=0;
M_sys_con1(:,np_deg)=0;
np_z_diag=zeros(3*sNode_q,1);
np_z_diag(np_deg)=1;
np_z_diag_1=diag(np_z_diag);
M_sys_con1=M_sys_con1+np_z_diag_1;
% F_sys_con1(np_deg)=0;

function [K_sys,F_sys]=K_F_sys(sNode_q,Node,ncon,uabc,wabc,f0,q0,F_app,ax,ay,az)
global A E Ks G I sncon

K_sys=zeros(3*sNode_q,3*sNode_q);    F_sys=zeros(3*sNode_q,1);

for el=1:sncon
    n1=ncon(el,1);    n2=ncon(el,2);
    xa=Node(n1,1);    xc=Node(n2,1);
    
    wabc_e=wabc(el,:);     uabc_e=uabc(el,:);
    
    % K_ele=calculate_K_element(xa,xb,A,E,Ks,G,I,wab_e,ax,ay,az);
    % [K_ele T_ele]=calculate_K_T_element(xa,xc,A,E,Ks,G,I,wabc_e,uabc_e);
    [K_ele]=calculate_K__element(xa,xc,A,E,Ks,G,I,wabc_e,uabc_e,ax,ay,az);
    m0=0;
    F_ele=calculate_F_element(xa,xc,f0,q0,F_app,m0);
    % X_ele=x_ele(el,:);
    % R_ele=-(F_ele-K_ele*X_ele');
    
    % Making K_system and F_system
    r1=3*(n1-1+el-1)+1:3*(n1-1+el-1)+3;   % Note: Just for Unbranched frams
    r2=3*(n2-1+el-1)+1:3*(n2-1+el-1)+3;
    r3=3*(n2-1+el)+1:3*(n2-1+el)+3;
    
    rrr=[r1 r2 r3];
    % T_sys(rrr(:),rrr(:))=T_sys(rrr(:),rrr(:))+T_ele(:,:);
    % R_sys(rrr(:))=R_sys(rrr(:))+R_ele(:);
    
    K_sys(rrr(:),rrr(:))=K_sys(rrr(:),rrr(:))+K_ele(:,:);
    F_sys(rrr(:))=F_sys(rrr(:))+F_ele(:);
    
    
end


function [K_sys,F_sys]=K_F_sys1(sNode_q,Node,ncon,uabc,wabc,f1,q0,F_app,ax,ay,az)
global A E Ks G I sncon dx

K_sys=zeros(3*sNode_q,3*sNode_q);    F_sys=zeros(3*sNode_q,1);

for el=1:sncon
    n1=ncon(el,1);    n2=ncon(el,2);
    xa=Node(n1,1);    xc=Node(n2,1);
    
    wabc_e=wabc(el,:);     uabc_e=uabc(el,:);
    du=(uabc_e(1,1)-uabc_e(1,3))/dx;
    dw=(wabc_e(1,1)-wabc_e(1,3))/dx;
    dww=(wabc_e(1,1)-wabc_e(1,3))/dx^2;
    N=E*A*(du+(0.5*dw^2));
    m0=E*I*(dww/(1+dw^2));
    f0=N+f1;

    
    % K_ele=calculate_K_element(xa,xb,A,E,Ks,G,I,wab_e,ax,ay,az);
    % [K_ele T_ele]=calculate_K_T_element(xa,xc,A,E,Ks,G,I,wabc_e,uabc_e);
    [K_ele]=calculate_K__element(xa,xc,A,E,Ks,G,I,wabc_e,uabc_e,ax,ay,az);
    
    F_ele=calculate_F_element(xa,xc,f0,q0,F_app,m0);
    % X_ele=x_ele(el,:);
    % R_ele=-(F_ele-K_ele*X_ele');
    
    % Making K_system and F_system
    r1=3*(n1-1+el-1)+1:3*(n1-1+el-1)+3;   % Note: Just for Unbranched frams
    r2=3*(n2-1+el-1)+1:3*(n2-1+el-1)+3;
    r3=3*(n2-1+el)+1:3*(n2-1+el)+3;
    
    rrr=[r1 r2 r3];
    % T_sys(rrr(:),rrr(:))=T_sys(rrr(:),rrr(:))+T_ele(:,:);
    % R_sys(rrr(:))=R_sys(rrr(:))+R_ele(:);
    
    K_sys(rrr(:),rrr(:))=K_sys(rrr(:),rrr(:))+K_ele(:,:);
    F_sys(rrr(:))=F_sys(rrr(:))+F_ele(:);
    
    
end

function [K_sys,F_sys]=K_F_sys2(sNode_q,Node,ncon,uabc,wabc,f1,q0,F_app,ax,ay,az,F)
global A E Ks G I sncon dx

K_sys=zeros(3*sNode_q,3*sNode_q);    F_sys=zeros(3*sNode_q,1);
i=1;
j=2;
for el=1:sncon
    n1=ncon(el,1);    n2=ncon(el,2);
    xa=Node(n1,1);    xc=Node(n2,1);
    
    wabc_e=wabc(el,:);     uabc_e=uabc(el,:);

    f=F(i,1);
    m0=F(j,1);
    f0=f1+f;
    i=i+6;
    j=j+6;
    
    % K_ele=calculate_K_element(xa,xb,A,E,Ks,G,I,wab_e,ax,ay,az);
    % [K_ele T_ele]=calculate_K_T_element(xa,xc,A,E,Ks,G,I,wabc_e,uabc_e);
    [K_ele]=calculate_K__element(xa,xc,A,E,Ks,G,I,wabc_e,uabc_e,ax,ay,az);
    
    F_ele=calculate_F_element(xa,xc,f0,q0,F_app,m0);
    % X_ele=x_ele(el,:);
    % R_ele=-(F_ele-K_ele*X_ele');
    
    % Making K_system and F_system
    r1=3*(n1-1+el-1)+1:3*(n1-1+el-1)+3;   % Note: Just for Unbranched frams
    r2=3*(n2-1+el-1)+1:3*(n2-1+el-1)+3;
    r3=3*(n2-1+el)+1:3*(n2-1+el)+3;
    
    rrr=[r1 r2 r3];
    % T_sys(rrr(:),rrr(:))=T_sys(rrr(:),rrr(:))+T_ele(:,:);
    % R_sys(rrr(:))=R_sys(rrr(:))+R_ele(:);
    
    K_sys(rrr(:),rrr(:))=K_sys(rrr(:),rrr(:))+K_ele(:,:);
    F_sys(rrr(:))=F_sys(rrr(:))+F_ele(:);
    
    
end


function [M_sys]=Mass_sys1(sNode_q,Node,ncon,am)
global A  I sncon ro_t

M_sys=zeros(3*sNode_q,3*sNode_q);
% F_sys=zeros(3*sNode_q,1);

for el=1:sncon
    n1=ncon(el,1);    n2=ncon(el,2);
    xa=Node(n1,1);    xc=Node(n2,1);
    
    %     wabc_e=wabc(el,:);     uabc_e=uabc(el,:);
    
    % K_ele=calculate_K_element(xa,xb,A,E,Ks,G,I,wab_e,ax,ay,az);
    % [K_ele T_ele]=calculate_K_T_element(xa,xc,A,E,Ks,G,I,wabc_e,uabc_e);
    %     [K_ele]=calculate_K__element(xa,xc,A,E,Ks,G,I,wabc_e,uabc_e,ax,ay,az);
    
    %     F_ele=calculate_F_element(xa,xc,f0,q0,F_app);
    % X_ele=x_ele(el,:);
    % R_ele=-(F_ele-K_ele*X_ele');
    
    % Making K_system and F_system
    r1=3*(n1-1+el-1)+1:3*(n1-1+el-1)+3;   % Note: Just for Unbranched frams
    r2=3*(n2-1+el-1)+1:3*(n2-1+el-1)+3;
    r3=3*(n2-1+el)+1:3*(n2-1+el)+3;
    
    rrr=[r1 r2 r3];
    %     % T_sys(rrr(:),rrr(:))=T_sys(rrr(:),rrr(:))+T_ele(:,:);
    %     % R_sys(rrr(:))=R_sys(rrr(:))+R_ele(:);
    %
    %     K_sys(rrr(:),rrr(:))=K_sys(rrr(:),rrr(:))+K_ele(:,:);
    %     F_sys(rrr(:))=F_sys(rrr(:))+F_ele(:);
    
    
    
    [M_ele]=calculate_M_element(xa,xc,ro_t,A,I);
    
    M_sys(rrr(:),rrr(:))=M_sys(rrr(:),rrr(:))+(am*M_ele(:,:));
    
    
end



function [M_sys]=Mass_sys2(sNode_q,Node,ncon,am)
global A  I sncon ro_t

M_sys=zeros(3*sNode_q,3*sNode_q);
% F_sys=zeros(3*sNode_q,1);

for el=1:sncon
    n1=ncon(el,1);    n2=ncon(el,2);
    xa=Node(n1,1);    xc=Node(n2,1);
    
    %     wabc_e=wabc(el,:);     uabc_e=uabc(el,:);
    
    % K_ele=calculate_K_element(xa,xb,A,E,Ks,G,I,wab_e,ax,ay,az);
    % [K_ele T_ele]=calculate_K_T_element(xa,xc,A,E,Ks,G,I,wabc_e,uabc_e);
    %     [K_ele]=calculate_K__element(xa,xc,A,E,Ks,G,I,wabc_e,uabc_e,ax,ay,az);
    
    %     F_ele=calculate_F_element(xa,xc,f0,q0,F_app);
    % X_ele=x_ele(el,:);
    % R_ele=-(F_ele-K_ele*X_ele');
    
    % Making K_system and F_system
    r1=3*(n1-1+el-1)+1:3*(n1-1+el-1)+3;   % Note: Just for Unbranched frams
    r2=3*(n2-1+el-1)+1:3*(n2-1+el-1)+3;
    r3=3*(n2-1+el)+1:3*(n2-1+el)+3;
    
    rrr=[r1 r2 r3];
    %     % T_sys(rrr(:),rrr(:))=T_sys(rrr(:),rrr(:))+T_ele(:,:);
    %     % R_sys(rrr(:))=R_sys(rrr(:))+R_ele(:);
    %
    %     K_sys(rrr(:),rrr(:))=K_sys(rrr(:),rrr(:))+K_ele(:,:);
    %     F_sys(rrr(:))=F_sys(rrr(:))+F_ele(:);
    
    
    
    [M_ele]=mass(xa,xc,ro_t,A,I);
    
    M_ele=am*M_ele;
    
    M_sys(rrr(:),rrr(:))=M_sys(rrr(:),rrr(:))+(am*M_ele(:,:));
    
    
end


function [M_sys]=Mass_sys3(sNode_q,Node,ncon,am,dx)
global A  I sncon ro_t

M_sys=zeros(3*sNode_q,3*sNode_q);
% F_sys=zeros(3*sNode_q,1);

for el=1:sncon
    n1=ncon(el,1);    n2=ncon(el,2);
    xa=Node(n1,1);    xc=Node(n2,1);
    
    %     wabc_e=wabc(el,:);     uabc_e=uabc(el,:);
    
    % K_ele=calculate_K_element(xa,xb,A,E,Ks,G,I,wab_e,ax,ay,az);
    % [K_ele T_ele]=calculate_K_T_element(xa,xc,A,E,Ks,G,I,wabc_e,uabc_e);
    %     [K_ele]=calculate_K__element(xa,xc,A,E,Ks,G,I,wabc_e,uabc_e,ax,ay,az);
    
    %     F_ele=calculate_F_element(xa,xc,f0,q0,F_app);
    % X_ele=x_ele(el,:);
    % R_ele=-(F_ele-K_ele*X_ele');
    
    % Making K_system and F_system
    r1=3*(n1-1+el-1)+1:3*(n1-1+el-1)+3;   % Note: Just for Unbranched frams
    r2=3*(n2-1+el-1)+1:3*(n2-1+el-1)+3;
    r3=3*(n2-1+el)+1:3*(n2-1+el)+3;
    
    rrr=[r1 r2 r3];
    %     % T_sys(rrr(:),rrr(:))=T_sys(rrr(:),rrr(:))+T_ele(:,:);
    %     % R_sys(rrr(:))=R_sys(rrr(:))+R_ele(:);
    %
    %     K_sys(rrr(:),rrr(:))=K_sys(rrr(:),rrr(:))+K_ele(:,:);
    %     F_sys(rrr(:))=F_sys(rrr(:))+F_ele(:);
    
    
    
    [M_ele]=mass_mms(xa,xc,ro_t,A,I,dx);
    
    
    M_sys(rrr(:),rrr(:))=M_sys(rrr(:),rrr(:))+(am*M_ele(:,:));
    
    
end


function [K_ele]=calculate_K__element(xa,xc,A,E,Ks,G,I,wabc_e,uabc_e,ax,ay,az)
Axx=E*A;
Dxx=E*I;
Sxx=Ks*G*A;
xb=0.5*(xa+xc);
% syms x L F a b n  MT R  c s  xi real ;    %   b > a
% syms xa xb
% xab=[xa xb];
% S1(1)=0.5*(1-xi);  S1(2)=0.5*(1+xi)
% S2(1)=0.5*(1-xi);  S2(2)=0.5*(1+xi);
% S3(1)=0.5*(1-xi);  S3(2)=0.5*(1+xi);

% x=xab(1)*S1(1)+xab(2)*S1(2);


% Dx_xi=diff(x,xi);   Jac=Dx_xi;
% for i=1:2
% DSx(i)=diff(S1(i),xi)/Jac;
% DSx(i)=simplify(DSx(i));
% end

i=1:3 ; j=1:3;
% M=[];
% M=int((S1(i))'*S1(j),xi,-1,1)*Jac;
% M=simplify(M);
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
% F_T12_extra=0;    % the word extra is introduced here because T12 has other terms: K12
% F_T22_extra_1=0;
% F_T22_extra_2=0;
for ite=1:LGP
    xi=Gauss_reduced(ite,1);  wt=Gauss_reduced(ite,2);
    S2(1)=-0.5*(1-xi)*xi;  S2(2)=(1-xi*xi);    S2(3)=0.5*(1+xi)*xi;   %S1 =[ 1/2-1/2*xi, 1/2+1/2*xi];
    Jac=xa*xi-1/2*xa-2*xb*xi+xc*xi+1/2*xc;   % Jac =-1/2*xa+1/2*xb;
    DSx =[ -(2*xi-1)/(-2*xa*xi+xa+4*xb*xi-2*xc*xi-xc),      4*xi/(-2*xa*xi+xa+4*xb*xi-2*xc*xi-xc), -(2*xi+1)/(-2*xa*xi+xa+4*xb*xi-2*xc*xi-xc)];
    %DSx =[1/(xa-xb), -1/(xa-xb)];
    %     f_T12_extra=(wa*DSx(1)+wb*DSx(2)).*((DSx(i))'*DSx(j));
    %     F_T12_extra=F_T12_extra+f_T12_extra*wt*Jac;
    
    %     f_T22_extra_1=(ua*DSx(1)+ub*DSx(2)+uc*DSx(3)).*((DSx(i))'*DSx(j));
    %     F_T22_extra_1=F_T22_extra_1+f_T22_extra_1*wt*Jac;
    
    %     f_T22_extra_2=((wa*DSx(1)+wb*DSx(2)+wc*DSx(3))^2).*(DSx(i))'*DSx(j);
    %     F_T22_extra_2=F_T22_extra_2+f_T22_extra_2*wt*Jac;
    
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

%  T11=K11;
%  T12=K12+0.5*Axx*F_T12_extra;    %T12=2*K12;
%  T12=2*K12;
%  T21=K21;
%  T13=K13;
%  T31=K31;
%  T22=K22+Axx*(F_T22_extra_1+F_T22_extra_2);
%  T23=K23;
%  T32=K32;
%  T33=K33;

%  T_ele([1 4 7],[1 4 7])=T11;
%  T_ele([2 5 8],[2 5 8])=T22;
%  T_ele([3 6 9],[3 6 9])=T33;
%  T_ele([2 5 8],[1 4 7])=T21;
%  T_ele([1 4 7],[2 5 8])=T12;
%  T_ele([1 4 7],[3 6 9])=T13;
%  T_ele([2 5 8],[3 6 9])=T23;
%  T_ele([3 6 9],[1 4 7])=T31;
%  T_ele([3 6 9],[2 5 8])=T32;

return

%% function calculate_K_T_element
function [K_ele, T_ele]=calculate_K_T_element(xa,xc,A,E,Ks,G,I,wabc_e,uabc_e,ax,ay,az)%non
Axx=E*A;
Dxx=E*I;
Sxx=Ks*G*A;
xb=0.5*(xa+xc);
% ax=1;   ay=1;   az=1;
% syms x L F a b n  MT R  c s  xi real ;    %   b > a
% syms xa xb
% xab=[xa xb];
% S1(1)=0.5*(1-xi);  S1(2)=0.5*(1+xi)
% S2(1)=0.5*(1-xi);  S2(2)=0.5*(1+xi);
% S3(1)=0.5*(1-xi);  S3(2)=0.5*(1+xi);

% x=xab(1)*S1(1)+xab(2)*S1(2);


% Dx_xi=diff(x,xi);   Jac=Dx_xi;
% for i=1:2
% DSx(i)=diff(S1(i),xi)/Jac;
% DSx(i)=simplify(DSx(i));
% end

i=1:3 ; j=1:3;
% M=[];
% M=int((S1(i))'*S1(j),xi,-1,1)*Jac;
% M=simplify(M);
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
    
    f_T22_extra_1=ax*ay*(ua*DSx(1)+ub*DSx(2)+uc*DSx(3)).*((DSx(i))'*DSx(j));
    F_T22_extra_1=F_T22_extra_1+f_T22_extra_1*wt*Jac;
    
    f_T22_extra_2=ay*((wa*DSx(1)+wb*DSx(2)+wc*DSx(3))^2).*(DSx(i))'*DSx(j);
    F_T22_extra_2=F_T22_extra_2+f_T22_extra_2*wt*Jac;
    
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

T11=K11;
%  T12=K12+0.5*Axx*F_T12_extra;    %T12=2*K12;
T12=2*K12;
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

return

%% function calculate_K_element
function K_ele=calculate_K_element(xa,xb,A,E,Ks,G,I,wab_e,ax,ay,az)
Axx=E*A;
Dxx=E*I;
Sxx=Ks*G*A;
% syms x L F a b n  MT R  c s  xi real ;    %   b > a
% syms xa xb
xab=[xa xb];
% S1(1)=0.5*(1-xi);  S1(2)=0.5*(1+xi)
% S2(1)=0.5*(1-xi);  S2(2)=0.5*(1+xi);
% S3(1)=0.5*(1-xi);  S3(2)=0.5*(1+xi);

% x=xab(1)*S1(1)+xab(2)*S1(2);
% ax=1;   ay=1;   az=1;

% Dx_xi=diff(x,xi);   Jac=Dx_xi;
% for i=1:2
% DSx(i)=diff(S1(i),xi)/Jac;
% DSx(i)=simplify(DSx(i));
% end

i=1:2 ; j=1:2;
% M=[];
% M=int((S1(i))'*S1(j),xi,-1,1)*Jac;
% M=simplify(M);
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
    
    f_K11=ax*(DSx(i))'*DSx(j);
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
    
    f_K12=ax*ay*(wa*DSx(1)+wb*DSx(2)).*((DSx(i))'*DSx(j));
    F_K12=F_K12+f_K12*wt*Jac;
    
    f_K22_1=az*(DSx(i))'*DSx(j);
    F_K22_1=F_K22_1+f_K22_1*wt*Jac;
    
    f_K22_2=ax*ay*((wa*DSx(1)+wb*DSx(2))^2).*(DSx(i))'*DSx(j);
    F_K22_2=F_K22_2+f_K22_2*wt*Jac;
    
    f_K23=az*(DSx(i))'*S1(j);   % or maybe this is true:   DSx(i)*(S1(j))';
    F_K23=F_K23+f_K23*wt*Jac;
    
    f_K33_1=ay*(DSx(i))'*DSx(j);
    F_K33_1=F_K33_1+f_K33_1*wt*Jac;
    
    f_K33_2=az*(S1(i))'*S1(j);
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

%% function  calculate_F_element
function F_ele=calculate_F_element(xa,xc,f0,q0,F0,m0)
% F_ele=zeros(9,1);




% syms xi
% S1(1)=0.5*(1-xi);  S1(2)=0.5*(1+xi);
% S2(1)=0.5*(1-xi);  S2(2)=0.5*(1+xi);
% S3(1)=0.5*(1-xi);  S3(2)=0.5*(1+xi);

% x=xab(1)*S1(1)+xab(2)*S1(2);


% Dx_xi=diff(x,xi);   Jac=Dx_xi;
% % for i=1:2
% % DSx(i)=diff(S1(i),xi)/Jac;
% % DSx(i)=simplify(DSx(i));
% % end
xb=0.5*(xa+xc);
i=1:3 ; j=1:3;
% M=[];
% M=int((S1(i))'*S1(j),xi,-1,1)*Jac;
% M=simplify(M);
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
% F3=zeros(3,1);    
F3=m0*F_F2;
F_ele([1;4;7])=F1;
F_ele([2;5;8])=F2;
F_ele([3;6;9])=F3;
F_ele=F_ele';     %F_ele=double(F_ele');
% F_ele=vpa(F_ele');
if xa==0
    F_ele(1,1)=F0;
end
return








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
consTy=sqrt(4.8);
gauss{4}=[sqrt((3+consTy)/7) 0.5-1/(3*consTy)
    sqrt((3-consTy)/7) 0.5+1/(3*consTy)
    -sqrt((3-consTy)/7) 0.5+1/(3*consTy)
    -sqrt((3+consTy)/7) 0.5-1/(3*consTy)];
gauss{5}=[0.5384693101 0.4786286705
    0.9061798459 0.2369268850
    0.0000000000 0.5688888889
    -0.5384693101 0.4786286705
    -0.9061798459 0.2369268850];
Gauss=gauss{NGP};
return

%% function input4
function [Node ncon Constraints]=input4()

c1=1e27;
c2=1e-27;
%[node_number rigidity_u rigidity_w rigidity_teta]
Constraints=[1	c2	c1	c2
    26	c1	c2	c1];

force=[];

Node=[
    0	0	0
    2	0	0
    4	0	0
    6	0	0
    8	0	0
    10	0	0
    12	0	0
    14	0	0
    16	0	0
    18	0	0
    20	0	0
    22	0	0
    24	0	0
    26	0	0
    28	0	0
    30	0	0
    32	0	0
    34	0	0
    36	0	0
    38	0	0
    40	0	0
    42	0	0
    44	0	0
    46	0	0
    48	0	0
    50	0	0
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
    21	22	0.797	0.00004	30000000	12000000	7900	0.008	0.3
    22	23	0.797	0.00004	30000000	12000000	7900	0.008	0.3
    23	24	0.797	0.00004	30000000	12000000	7900	0.008	0.3
    24	25	0.797	0.00004	30000000	12000000	7900	0.008	0.3
    25	26	0.797	0.00004	30000000	12000000	7900	0.008	0.3
    ];

return

%% function input3
function [Node ncon Constraints]=input3()

c1=1e27;
c2=1e-27;
%[node_number rigidity_u rigidity_w rigidity_teta]
Constraints=[1	c2	c1	c2
    9	c1	c2	c1];

force=[];

%3: inner diameter 4:epsilon=raphness   5:E  6:G   7:Rho   8:thickness  9:poisson ratio

Node=[
    0	0	0
    6.25	0	0
    12.5	0	0
    18.75	0	0
    25	0	0
    31.25	0	0
    37.5	0	0
    43.75	0	0
    50	0	0
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
    ];

return




%% function input2
function [Node ncon Constraints]=input2()

c1=1e270;
c2=1e-270;
%[node_number rigidity_u rigidity_w rigidity_teta]
Constraints=[1	c2	c1	c2
    17	c1	c2	c1];

force=[];

%3: inner diameter 4:epsilon=raphness   5:E  6:G   7:Rho   8:thickness  9:poisson ratio
Node=[
    0	0	0
    3.125	0	0
    6.25	0	0
    9.375	0	0
    12.5	0	0
    15.625	0	0
    18.75	0	0
    21.875	0	0
    25	0	0
    28.125	0	0
    31.25	0	0
    34.375	0	0
    37.5	0	0
    40.625	0	0
    43.75	0	0
    46.875	0	0
    50	0	0
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
    ];


return

%% function inpuTy
function [Node ncon Constraints]=inpuTy()
% global Reaction1 Reaction2 X_ma Xd_ma Xdd_ma  x_loc3d  xd_loc3d  xdd_loc3d
% global Coupled_Analys Ttotal
% global L_deltax      c_glo      for_newmark    no_or_with_damping


% no_or_with_damping=1 ;%structural damping, associate with matrix C is the aim of this part
% omega_1=40;

% c_glo=1209.55;       c=c_glo;
% deltax= 2;        Ttotal=2;      dt=deltax/c;
% for_newmark=1;   % such as rigidity, this parameter can be discussed for diffrent values
% max_nwhile=40;   %i suggest that use this to be as 100 but write the appropriate nomber in a matrix namely "iterate_each_step" it concerns me that it falls in an infinite loop. if the mentioned matrix has a value more that 80, this guess must be checked for that time step


% nomber_of_deformed_shapes=1;
% elti=Ttotal/nomber_of_deformed_shapes;
% time=iiti:iiti:Ttotal-iiti;

% time=[0.05 0.1 0.15 0.2 0.25 0.3]



% PlotNod=[32t];

% col=' ';

%%Graph Options
%%%%PLOTPLOTPLOTPLOTPLOTPLOTPLOTPLOTPLOTPLOTPLOTPLOTPLOTPLOTPLOTPLOTPLOT


% new_figure=0;%1
% PlotNod_h=[166];
% % PlotNod_Q=[1 17 21];
% Reaction_Node=[1 166];
% X_glo_Node=[156];
% % X_loc_Ele_nomber=[8 16 20];  %the nomber of elements is equal to the nomber of nodes minus one
% % Force_loc_Ele_nomber=[8 16 20]; %Internal Forces

%%%%PLOTPLOTPLOTPLOTPLOTPLOTPLOTPLOTPLOTPLOTPLOTPLOTPLOTPLOTPLOTPLOTPLOT
% h   plot   :  1
% Q  plot   :  2
% Plot_Setting=0;
%%Friction Options
% withORwithoutFriction=0;  % ANALYSIS WITH OR WITHOUT DAMPING(FRICTION)

%%zero
% the thirth column is rigidity, if a point is quite free to move, you may
% write a very low value for it or you can do not write any nomber of that
% point
% const_1=1e200;
% const_2=1e-200;
% zero=[
%     1 1 const_1
%     1 2 const_1
%     1 3 const_1
%     156 1 const_2
%     156 2 const_2
%     156 3 const_2
%     166 1 const_1
%     166 2 const_1
%     166 3 const_1
%     ];
c1=1e270;
c2=1e-270;
%[node_number rigidity_u rigidity_w rigidity_teta]
Constraints=[1	c2	c1	c2
    5	c1	c2	c1];

force=[];

%%Node  &   ncon
poisson_ratio=0 ; %No Pois. Coupl. poisson_ratio=0    Poi. Coup. included    =0.3

%%ncon
%3: inner diameter 4:epsilon=raphness   5:E  6:G   7:Rho   8:thickness  9:poisson ratio
Node=[
    0	0	0
    12.5	0	0
    25	0	0
    37.5	0	0
    50	0	0
    ];
ncon=[
    1	2	0.797	0.00004	30000000	12000000	7900	0.008	0.3
    2	3	0.797	0.00004	30000000	12000000	7900	0.008	0.3
    3	4	0.797	0.00004	30000000	12000000	7900	0.008	0.3
    4	5	0.797	0.00004	30000000	12000000	7900	0.008	0.3
    ];


% ncon(:,9)=poisson_ratio;

% FV_mat=load('force_classical_1.txt');   %Dynamic force

% save('MatFileForStructuralAnalysis.mat')

return


%% function input5
function [Node ncon Constraints]=input5()

c1=1e270;
c2=1e-270;
%[node_number rigidity_u rigidity_w rigidity_teta]
% Constraints=[1	c2	c1	c2
% 5	c1	c2	c1];
Constraints=[1	c2	c1	c2
    9	c2	c1	c2];
force=[];

%%Node  &   ncon
poisson_ratio=0 ; %No Pois. Coupl. poisson_ratio=0    Poi. Coup. included    =0.3

%%ncon
%3: inner diameter 4:epsilon=raphness   5:E  6:G   7:Rho   8:thickness  9:poisson ratio
Node=[
    0	0	0
    12.5	0	0
    25	0	0
    37.5	0	0
    50	0	0
    62.5	0	0
    75	0	0
    87.5	0	0
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
    ];


% ncon(:,9)=poisson_ratio;

% FV_mat=load('force_classical_1.txt');   %Dynamic force

% save('MatFileForStructuralAnalysis.mat')

return

%% function input6
function [Node, ncon, Constraints]=input6()
global nsec


%[node_number rigidity_u rigidity_w rigidity_teta]
% Constraints=[1	c2	c1	c2
% 5	c1	c2	c1];

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
nsec=50;
% nsec=20;

Node=zeros(nsec+1,3);
Node(:,1)=0:6/nsec:6;
ncon=zeros(nsec,7);
ncon(:,1)=1:nsec;
ncon(:,2)=2:nsec+1;
c1=1e27;
c2=1e-27;
Constraints=[1	c2	c1	c1
    nsec+1	c1	c1	c1];
% Constraints=[1	c2	c1	c2
%     nsec+1	c1	c1	c2];
% ncon=[
%     1	2	0.797	0.00004	30000000	12000000	7900	0.008	0.3
%     2	3	0.797	0.00004	30000000	12000000	7900	0.008	0.3
%     3	4	0.797	0.00004	30000000	12000000	7900	0.008	0.3
%     4	5	0.797	0.00004	30000000	12000000	7900	0.008	0.3
%     5	6	0.797	0.00004	30000000	12000000	7900	0.008	0.3
%     6	7	0.797	0.00004	30000000	12000000	7900	0.008	0.3
%     7	8	0.797	0.00004	30000000	12000000	7900	0.008	0.3
%     8	9	0.797	0.00004	30000000	12000000	7900	0.008	0.3
%     9	10	0.797	0.00004	30000000	12000000	7900	0.008	0.3
%     10	11	0.797	0.00004	30000000	12000000	7900	0.008	0.3
%     11	12	0.797	0.00004	30000000	12000000	7900	0.008	0.3
%     12	13	0.797	0.00004	30000000	12000000	7900	0.008	0.3
%     13	14	0.797	0.00004	30000000	12000000	7900	0.008	0.3
%     14	15	0.797	0.00004	30000000	12000000	7900	0.008	0.3
%     15	16	0.797	0.00004	30000000	12000000	7900	0.008	0.3
%     16	17	0.797	0.00004	30000000	12000000	7900	0.008	0.3
%     17	18	0.797	0.00004	30000000	12000000	7900	0.008	0.3
%     18	19	0.797	0.00004	30000000	12000000	7900	0.008	0.3
%     19	20	0.797	0.00004	30000000	12000000	7900	0.008	0.3
%     20	21	0.797	0.00004	30000000	12000000	7900	0.008	0.3
%     ];


% ncon(:,9)=poisson_ratio;

% FV_mat=load('force_classical_1.txt');   %Dynamic force

% save('MatFileForStructuralAnalysis.mat')

return












syms x

F(1)=1-3*(x/L)^2+2*(x/L)^3;   % "F" is substituted for  Fi
F(2)=-x*(1-x/L)^2;
F(3)=3*(x/L)^2-2*(x/L)^3;
F(4)=-x*((x/L)^2-x/L);
i=[1:4] ; j=[1:4];
M2=int((F(i))'*F(j),x,0,L);
% MM=(M2.*420)./L;
%  C1=int((diff(F(i),x))'*diff(F(j),x),x,0,L);
%  K1=int((diff(diff(F(i),x),x))'*diff(diff(F(j),x),x),x,0,L);
a=0;b=L;
N(1)=((b-x)/(b-a));
N(2)=((a-x)/(a-b));
ii=[1:2] ; jj=[1:2];
M1=int((N(ii))'*N(jj),x,a,b);

iM1=[1 4];
MT(1:6,1:6)=0 ;  %=zeros(6,6)      %  MT : Mass Total
MT(iM1,iM1)=M1;

iM2=[2 3 5 6];
MT(iM2,iM2)=M2;

MTfac=MT./(L/420);
simple(MTfac);

R(1:6,1:6)=0;
R(3,3)=1;  R(6,6)=1;
R(1,1)=c;R(2,2)=c;R(1,2)=s;R(2,1)=-s;
R(4,4)=c;R(5,5)=c;R(4,5)=s;R(5,4)=-s;

MTG=R'*MT*R;  %  MTG : Mass Total Global
% c=1;s=0;
% mmat=eval(MTG)
MTGfac=MTG./(L/420);
MTGfac=simple(MTGfac);
%  s11=simplify(s1); s1=expand(s1) ; s1=simplify(s1); s1=factor(s1) ; s11=horner(s11); pretty(s11)
% pause;
clear;
syms ke E I L c s R G J  A   cc  real;
ke(1:6,1:6)=0 ;  %=zeros(6,6)

RR=A*L^2/I    %in plane element I=Iz
ke(1,1)=RR;
%  ke(1,2)=0;
ke(2,2)=12;
% ke(1,3)=0;
ke(2,3)=6*L;
ke(3,3)=4*L^2;
ke(3,6)=2*L^2;
%----------------------------------------------
ke(1,4)=-ke(1,1);
%   ke(1,5)=-ke(1,2);        because the element is placed on "ox" cordinate
%   ke(1,6)=ke(1,3);
%   ke(2,4)=-ke(1,2);
ke(2,5)=-ke(2,2);
ke(2,6)=ke(2,3);
%   ke(3,4)=-ke(1,3);
ke(3,5)=-ke(2,3);
ke(4,4:5)=ke(1,1:2);
%   ke(4,6)=-ke(1,3);
ke(5,5)=ke(2,2);
ke(5,6)=-ke(2,6);
ke(6,6)=ke(3,3);
%   ke=ke*E*I/L^3;
ke=ke+ke'-diag(diag(ke));
%   ke=ke*E*I/L^3;
R(1:6,1:6)=0;
R(3,3)=1;  R(6,6)=1;
R(1,1)=c;R(2,2)=c;R(1,2)=s;R(2,1)=-s;
R(4,4)=c;R(5,5)=c;R(4,5)=s;R(5,4)=-s;

k=R'*ke*R;
k=simple(k)
% kf=k./(E/L)
%  kf=k./(E*I/L^3)
kf=simplify(kf)
kf_factor=factor(kf)
pretty(kf)

%% function Mass
function [m]=mass(xa,xc,ro,A,I);
    
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
 
 S3=[1-(23*x^2)+(66*x^3)-(68*x^4)+(24*x^5)
    x*(1-(6*x)+(13*x^2)-(12*x^3)+(4*x^4))
    (16*x^2)*(1-(2*x)+(x^2))
    (-8*x^2)*(1-(4*x)+(5*x^2)-(2*x^3))
    (x^2)*(7-(34*x)+(52*x^2)-(24*x^3))
    (-x^2)*(1-(5*x)+(8*x^2)-(4*x^3))];

    syms x
for i=1:6
    for j=1:6
    m22(i,j)=int(S3(i)*S3(j),0,1);
    end
end
 m11=m11*ro*A;
 m22=m22*ro*I;
 m_ele([1 4 7],[1 4 7])=m11;
 m_ele([2 5 8],[2 5 8])=m22([1 3 5],[1 3 5]);
 m_ele([3 6 9],[3 6 9])=m22([2 4 6],[2 4 6]);
 m_ele([2 5 8],[1 4 7])=0;%m21;
 m_ele([1 4 7],[2 5 8])=0;%m12;
 m_ele([1 4 7],[3 6 9])=0;%m13;
 m_ele([2 5 8],[3 6 9])=m22([1 3 5],[2 4 6]);
 m_ele([3 6 9],[1 4 7])=0;%m31;
 m_ele([3 6 9],[2 5 8])=m22([2 4 6],[1 3 5]);
 m= m_ele;
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

function [M_ele]=calculate_M_element(xa,xc,ro_t,A,I)
% Axx=E*A;
% Dxx=E*I;
% Sxx=Ks*G*A;
RoA=ro_t*A;
RoI=ro_t*I;
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

function [X,Xd,Xdd] =Newmark(M,C,K,F,dt,Xi,Xdi,Xddi,sNode_q,n,nty)

persistent gamma_ beta_
if isempty(beta_)
    gamma_=1/2;              beta_=1/4;
end
c0 = 1/(beta_*dt*dt) ;
c1 = gamma_/(beta_*dt) ;
c2 = 1/(beta_*dt) ;
c3 = 1/(beta_*2) - 1 ;
c4 = gamma_/beta_ - 1 ;
% c5 = 0.5*dt*(gamma_/beta_ - 2 ) ;
c5 = dt*(gamma_/beta_ - 1 ) ;
c6 = dt*(1 - gamma_ ) ;
c7 = dt* gamma_  ;
Keff = c0*M + c1*C + K ;
Feff = F+ M*(c0*Xi+c2*Xdi+c3*Xddi) ...
    +C*(c1*Xi+c4*Xdi+c5*Xddi) ;

% if n<nty
%     n=n*3+1;
%     
%     g=n:1:3*sNode_q;
%     g=g';
%     Keff(g,:)=0;
%     Keff(:,g)=0;
%     np_z_diag=zeros(3*sNode_q,1);
%     np_z_diag(g)=1;
%     np_z_diag_1=diag(np_z_diag);
%     Keff=Keff+np_z_diag_1;
%     Feff(g)=0;
% end

% warning off all

% X=tridiagonal_matrix_algorithm(Keff,Feff);

X=Keff\Feff;
% X=X+Xi;

% [L_meth,U_meth] = lu(sparse(Keff));
% X = U_meth\(L_meth\Feff);

% X=pcg(sparse(Keff), Feff);


% fid = fopen('Xfile.txt', 'wt');
% fprintf(fid,'%3d\n', X);
% fclose(fid)

Xdd= c0*(X-Xi) - c2*Xdi - c3*Xddi ;
Xd = Xdi + c6*Xddi + c7*Xdd ;
% warning on all



function X=tridiagonal_matrix_algorithm(A,d)
persistent sn
if isempty(sn)%~exist('sn')
    sn=size(A,1);
end
a=diag(A,-1); a=[0;a];   b=diag(A);  c=diag(A,1);

for k = 2:sn
    m = a(k)/b(k - 1);
    b(k) = b(k) - m*c(k - 1);
    d(k) = d(k) - m*d(k - 1);
end

% Backward substitution phase
X(sn) = d(sn)/b(sn);
for k = sn-1:-1:1
    X(k) = (d(k) - c(k)* X(k + 1)) /b(k);
end
X=X';


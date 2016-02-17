%test1 Nonlinear dynamic analysis of Timoshenko beams by BEM.
%Part II: applications and validation
%Shear deformation effect in nonlinear free vibration analysis
%A = 150 × 10?6 m2          ? = 2.7 tn/m3
% E = 70 GPa                        G = 27 GPa
% Iy = 7812.5 × 10?12 m4 Iz = 450 × 10?12 m4
% ay = 1.74                             az = 1.20
% investigation for axial loads
%    The Plane Frame Element K M  Global Matrises

function []=static_nonlinear()
clc
clear var all
[Node, ncon, Constraints]= input8();
%% Preliminaries of the main while loop

global A E Ks G I sncon ro_t nsec initial_condition ...
    Xmat Xdmat Xddmat toler iter max_iter np_deg coef_damp ...
    n M_sys_con1 dt Load t nt dx A Ay Az I Pcr nx  ...
     TT nu  L st_analysis_type r

% Axx=E*A;
% Dxx=E*I;
% Sxx=Ks*G*A;
% xa=0;    xb=12.5;%syms xa xb




sncon=size(ncon,1);
sNode=size(Node,1);   sNode_q=2*sNode-1;
wabc=zeros(sncon,3);      uabc=zeros(sncon,3);
% K_sys=zeros(3*sNode,3*sNode);    F_sys=zeros(3*sNode,1);



X_sys_former=zeros(3*sNode_q,1);
x_ele=zeros(sncon,9);  %10*ones(sncon,6);
% % Applying Boundary Conditions, Method1
% np= repmat(3*(Constraints(:,1)-1),1,3)+repmat(1:3,size(Constraints,1),1);
% np=reshape(np',[],1);    %np=unique(a1);
% % making the rigidity vector
% Const_rigid=Constraints;   Const_rigid(:,1)=[];
% rigidity=zeros(3*sNode,1);
% rigidity(np)=reshape(Const_rigid',[],1);
% K_1=diag(rigidity);
% -----------
%  Applying Boundary Conditions, Method2
[as bs]=find(Constraints>1000);
np_deg=sort(2*3*(-1+Constraints(as))+bs-ones(size(bs,1),1));
% ----------------------
% wab_former=zeros(sncon,1);
%% LOADING
f0=0;
% q=0.00000000000000001;
q=1;
% load_increment_q0=0.00000000000000001;
load_increment_q0=0;
% BUCKLING=2466;

% F0_=BUCKLING;

% load_increment_F0=1;


% NI_F0=F0_/load_increment_F0;
% inp_graph=zeros(NI_F0,1);
% F0=load_increment_F0;
% Determinant=zeros(NI_F0,1);

% NI_q0=q/load_increment_q0;
q0=load_increment_q0;
% X_load_incre=zeros(3*sNode_q,NI_F0);
% F_app=0;

% persistent M_sys_con1
% M_sys=zeros(3*sNode_q,3*sNode_q);
Xmat=zeros(3*sNode_q,nt);
Xdmat=zeros(3*sNode_q,nt);
Xddmat=zeros(3*sNode_q,nt);
check=zeros(nt,2);


% if n==1&& iter==1
%     [M_ele]=calculate_M_element(xa,xc,ro_t,A,I);
%
%     M_sys(rrr(:),rrr(:))=M_sys(rrr(:),rrr(:))+M_ele(:,:);
% end

[M_sys]=Mass_sys(sNode_q,Node,ncon);
[M_sys_con1]=Mass_sys_con1(M_sys,np_deg,sNode_q);  % apply constraints Boundary conds)



results=zeros(nt,1);
    Load=zeros(1,nt);
%     Load(1,1)=0;
%     Load(1,2)=5;
%     Load(1,3)=8;
%     Load(1,4)=7;
%     Load(1,5)=5;
%     Load(1,6)=3;
%     Load(1,7)=2;
%     Load(1,8)=1;
%     Load(1,9)=0;
    


for n=1:nt-1
    %% The main while loop
    %         F_app=F_app+F0;   % F_app: applied load at the current load incerement
    F_app=Load(n+1)*1000;
%     F_app=100*sin(0.030*n*dt);
    
            [Fsys]=F_sys(sNode_q,Node,ncon,f0,q0,F_app);

    %     detlaF=F_app-Load(n);
    toler=100;   toler_w=100;   iter=0;
    X_sys=zeros(3*sNode_q,1);
    n
    wmax=n*3*r/nt;
            for j=1:(sNode_q);
                r1=3*(j-1)+1;    r2=3*(j-1)+2;      r3=3*(j-1)+3;
                j=j-1;
                X=dx*j/2;
                FI=pi;
                Xmat(r2,n)=wmax*sin(FI*X/L);
            end

    switch initial_condition
        case 'given'
            if n==1
                % % % %        for NE=20
                %                 X_sys=[0.0123370055063571,0,-1.81912751854593e-06,0.0120285803685739,2.74229132267169e-07,-1.81351686820713e-06,0.0117201552307937,5.46768409536662e-07,-1.79673097456298e-06,0.0114117300930196,8.15934942763568e-07,-1.76886194829885e-06,0.0111033049552544,1.08007352739942e-06,-1.73009286695601e-06,0.0107948798175006,1.33754967861163e-06,-1.68065172394364e-06,0.0104864546797607,1.58678359679647e-06,-1.62085410558548e-06,0.0101780295420366,1.82622944717116e-06,-1.55105827414327e-06,0.00986960440432997,2.05442173440242e-06,-1.47170455405620e-06,0.00956117926664194,2.26994133676877e-06,-1.38327265496421e-06,0.00925275412897325,2.47147314802453e-06,-1.28631679733496e-06,0.00894432899132413,2.65775970247156e-06,-1.18142632631987e-06,0.00863590385369433,2.82766866474593e-06,-1.06925571187217e-06,0.00832747871608313,2.98013518768664e-06,-9.50489423238684e-07,0.00801905357848937,3.11423758824740e-06,-8.25866064877575e-07,0.00771062844091144,3.22912985748167e-06,-6.96148376479440e-07,0.00740220330334735,3.32412365921534e-06,-5.62140910387302e-07,0.00709377816579475,3.39861265395135e-06,-4.24665895834903e-07,0.00678535302825104,3.45215880125420e-06,-2.84574022598859e-07,0.00647692789071338,3.48441036094272e-06,-1.42726763747623e-07,0.00616850275317877,3.49519037416721e-06,-4.90428180839087e-19,0.00586007761564416,3.48441036094287e-06,1.42726763746648e-07,0.00555165247810649,3.45215880125449e-06,2.84574022597903e-07,0.00524322734056278,3.39861265395179e-06,4.24665895834001e-07,0.00493480220301017,3.32412365921591e-06,5.62140910386468e-07,0.00462637706544607,3.22912985748235e-06,6.96148376478705e-07,0.00431795192786813,3.11423758824819e-06,8.25866064876954e-07,0.00400952679027435,2.98013518768752e-06,9.50489423238168e-07,0.00370110165266313,2.82766866474688e-06,1.06925571187178e-06,0.00339267651503331,2.65775970247256e-06,1.18142632631962e-06,0.00308425137738416,2.47147314802556e-06,1.28631679733485e-06,0.00277582623971543,2.26994133676979e-06,1.38327265496428e-06,0.00246740110202738,2.05442173440343e-06,1.47170455405645e-06,0.00215897596432072,1.82622944717211e-06,1.55105827414370e-06,0.00185055082659662,1.58678359679734e-06,1.62085410558608e-06,0.00154212568885668,1.33754967861239e-06,1.68065172394440e-06,0.00123370055110291,1.08007352740005e-06,1.73009286695690e-06,0.000925275413337634,8.15934942764052e-07,1.76886194829983e-06,0.000616850275563476,5.46768409536988e-07,1.79673097456402e-06,0.000308425137783268,2.74229132267334e-07,1.81351686820818e-06,0,0,1.81912751854699e-06]';
                %         X_sys=[0.0123215842451331,0,-6.54253206949244e-07,0.0120135446389887,9.86271122374053e-08,-6.52235312718465e-07,0.0117055050328447,1.96646460368480e-07,-6.46198198747426e-07,0.0113974654267015,2.93452787199434e-07,-6.36174998792709e-07,0.0110894258205594,3.88450793444134e-07,-6.22231562273573e-07,0.0107813862144188,4.81052632208096e-07,-6.04449892298841e-07,0.0104733466082800,5.70690125135695e-07,-5.82943494144915e-07,0.0101653070021433,6.56807309769577e-07,-5.57841221167501e-07,0.00985726739600876,7.38877118537779e-07,-5.29301439739378e-07,0.00954922778987668,8.16389162489752e-07,-4.97496680797579e-07,0.00924118818374711,8.88870462448466e-07,-4.62626272467749e-07,0.00893314857762007,9.55868769915520e-07,-4.24902175160768e-07,0.00862510897149553,1.01697683993973e-06,-3.84559770839631e-07,0.00831706936537339,1.07181169986893e-06,-3.41845230437794e-07,0.00800902975925351,1.12004186628136e-06,-2.97024194112195e-07,0.00770099015313568,1.16136307203504e-06,-2.50370982026370e-07,0.00739295054701964,1.19552775711877e-06,-2.02174951861888e-07,0.00708491094090509,1.22231785165472e-06,-1.52731818621438e-07,0.00677687133479168,1.24157581453675e-06,-1.02347532862242e-07,0.00646883172867906,1.25317514253103e-06,-5.13319232454024e-08,0.00616079212256682,1.25705219245931e-06,-3.13626497225405e-19,0.00585275251645458,1.25317514253112e-06,5.13319232447758e-08,0.00554471291034195,1.24157581453694e-06,1.02347532861633e-07,0.00523667330422853,1.22231785165500e-06,1.52731818620860e-07,0.00492863369811396,1.19552775711914e-06,2.02174951861361e-07,0.00462059409199789,1.16136307203548e-06,2.50370982025912e-07,0.00431255448588004,1.12004186628187e-06,2.97024194111817e-07,0.00400451487976014,1.07181169986949e-06,3.41845230437499e-07,0.00369647527363797,1.01697683994034e-06,3.84559770839434e-07,0.00338843566751340,9.55868769916144e-07,4.24902175160677e-07,0.00308039606138633,8.88870462449094e-07,4.62626272467762e-07,0.00277235645525672,8.16389162490370e-07,4.97496680797693e-07,0.00246431684912461,7.38877118538372e-07,5.29301439739588e-07,0.00215627724299010,6.56807309770132e-07,5.57841221167796e-07,0.00184823763685331,5.70690125136199e-07,5.82943494145292e-07,0.00154019803071448,4.81052632208537e-07,6.04449892299290e-07,0.00123215842457386,3.88450793444501e-07,6.22231562274085e-07,0.000924118818431748,2.93452787199719e-07,6.36174998793276e-07,0.000616079212288488,1.96646460368675e-07,6.46198198748039e-07,0.000308039606144442,9.86271122375047e-08,6.52235312719112e-07,0,0,6.54253206949899e-07]';
                X_sys=[0.0123191168438566,0,-5.58058957112973e-07,0.0120111389227485,8.41260588461734e-08,-5.56337749379084e-07,0.0117031610016407,1.67733711657152e-07,-5.51188259752993e-07,0.0113951830805334,2.50306691086836e-07,-5.42638751730702e-07,0.0110872051594270,3.31337223851364e-07,-5.30745394140863e-07,0.0107792272383217,4.10323893455084e-07,-5.15578134347433e-07,0.0104712493172176,4.86782060295111e-07,-4.97233789959617e-07,0.0101632713961151,5.60237504660785e-07,-4.75822269956474e-07,0.00985529347501417,6.30240653482009e-07,-4.51478656872787e-07,0.00954731555391502,6.96356160125086e-07,-4.24350115110618e-07,0.00923933763281770,7.58180587492156e-07,-3.94606666112199e-07,0.00893135971172222,8.15328181180942e-07,-3.62429106203399e-07,0.00862338179062856,8.67451573539965e-07,-3.28018209446210e-07,0.00831540386953665,9.14224100487129e-07,-2.91583952502266e-07,0.00800742594844638,9.55363016743096e-07,-2.53352918571286e-07,0.00769944802735760,9.90608791706057e-07,-2.13559094581255e-07,0.00739147010627012,1.01975026425514e-06,-1.72449293912216e-07,0.00708349218518373,1.04260142746069e-06,-1.30275752505693e-07,0.00677551426409817,1.05902790359333e-06,-8.72994368549741e-08,0.00646753634301318,1.06892178927926e-06,-4.37846213331603e-08,0.00615955842192847,1.07222879884423e-06,-2.11452137166759e-19,0.00585158050084377,1.06892178927933e-06,4.37846213327436e-08,0.00554360257975878,1.05902790359346e-06,8.72994368545732e-08,0.00523562465867321,1.04260142746088e-06,1.30275752505311e-07,0.00492764673758681,1.01975026425538e-06,1.72449293911860e-07,0.00461966881649932,9.90608791706351e-07,2.13559094580935e-07,0.00431169089541053,9.55363016743437e-07,2.53352918571011e-07,0.00400371297432025,9.14224100487509e-07,2.91583952502038e-07,0.00369573505322832,8.67451573540376e-07,3.28018209446044e-07,0.00338775713213464,8.15328181181374e-07,3.62429106203302e-07,0.00307977921103914,7.58180587492597e-07,3.94606666112174e-07,0.00277180128994179,6.96356160125524e-07,4.24350115110668e-07,0.00246382336884263,6.30240653482435e-07,4.51478656872913e-07,0.00215584544774171,5.60237504661187e-07,4.75822269956668e-07,0.00184786752663914,4.86782060295478e-07,4.97233789959877e-07,0.00153988960553509,4.10323893455407e-07,5.15578134347753e-07,0.00123191168442972,3.31337223851633e-07,5.30745394141235e-07,0.000923933763323277,2.50306691087044e-07,5.42638751731117e-07,0.000615955842215995,1.67733711657294e-07,5.51188259753438e-07,0.000307977921108142,8.41260588462456e-08,5.56337749379548e-07,0,0,5.58058957113443e-07]';
                % % %         for NE=40
                %         X_sys=[0.0123184999935466,0,-5.37094538384957e-07,0.0121645187436218,4.05141886188294e-08,-5.36680400095084e-07,0.0120105374936971,8.09659142763551e-08,-5.35438837149320e-07,0.0118565562437724,1.21292778704119e-07,-5.33371552960880e-07,0.0117025749938478,1.61432640354378e-07,-5.30481947805814e-07,0.0115485937439234,2.01323550087951e-07,-5.26774267766182e-07,0.0113946124939991,2.40904070493583e-07,-5.22254440073540e-07,0.0112406312440750,2.80113083478986e-07,-5.16929227301588e-07,0.0110866499941511,3.18890234532195e-07,-5.10807047095864e-07,0.0109326687442275,3.57175613402004e-07,-5.03897136830283e-07,0.0107786874943042,3.94910320562444e-07,-4.96210352669397e-07,0.0106247062443812,4.32036022951219e-07,-4.87758349838730e-07,0.0104707249944585,4.68495639145032e-07,-4.78554356093368e-07,0.0103167437445361,5.04232772894921e-07,-4.68612372655736e-07,0.0101627624946141,5.39192512788306e-07,-4.57947917322075e-07,0.0100087812446925,5.73320747073487e-07,-4.46577250993891e-07,0.00985479999477125,6.06565072878485e-07,-4.34518085841063e-07,0.00970081874485042,6.38873996548817e-07,-4.21788842195415e-07,0.00954683749493001,6.70197946817788e-07,-4.08409317404581e-07,0.00939285624501002,7.00488365587473e-07,-3.94399977668987e-07,0.00923887499509044,7.29698818804238e-07,-3.79782583463001e-07,0.00908489374517130,7.57783983288289e-07,-3.64579520681134e-07,0.00893091249525258,7.84700848465855e-07,-3.48814378770751e-07,0.00877693124533428,8.10407605493518e-07,-3.32511325310981e-07,0.00862294999541640,8.34864932438351e-07,-3.15695633292795e-07,0.00846896874549894,8.58034792545602e-07,-2.98393102986366e-07,0.00831498749558187,8.79881794943636e-07,-2.80630535117806e-07,0.00816100624566519,9.00371909463660e-07,-2.62435203589182e-07,0.00800702499574887,9.19473894480706e-07,-2.43835271634776e-07,0.00785304374583290,9.37157936208563e-07,-2.24859318644431e-07,0.00769906249591725,9.53397134874120e-07,-2.05536696733733e-07,0.00754508124600190,9.68166076876290e-07,-1.85897114587835e-07,0.00739109999608684,9.81442370131415e-07,-1.65970932247291e-07,0.00723711874617203,9.93205157898848e-07,-1.45788804537967e-07,0.00708313749625744,1.00343669383200e-06,-1.25381912255172e-07,0.00692915624634304,1.01212080663284e-06,-1.04781667377927e-07,0.00677517499642881,1.01924450509841e-06,-8.40198792261185e-08,0.00662119374651472,1.02479640306965e-06,-6.31285232764843e-08,0.00646721249660073,1.02876834457822e-06,-4.21398412684475e-08,0.00631323124668681,1.03115379879993e-06,-2.10861750470642e-08,0.00615924999677292,1.03194949528274e-06,-3.77855895653709e-19,0.00600526874685904,1.03115379879999e-06,2.10861750463110e-08,0.00585128749694513,1.02876834457833e-06,4.21398412676962e-08,0.00569730624703115,1.02479640306982e-06,6.31285232757329e-08,0.00554332499711707,1.01924450509863e-06,8.40198792253703e-08,0.00538934374720286,1.01212080663313e-06,1.04781667377184e-07,0.00523536249728849,1.00343669383234e-06,1.25381912254437e-07,0.00508138124737393,9.93205157899244e-07,1.45788804537244e-07,0.00492739999745914,9.81442370131865e-07,1.65970932246584e-07,0.00477341874754410,9.68166076876794e-07,1.85897114587156e-07,0.00461943749762878,9.53397134874675e-07,2.05536696733087e-07,0.00446545624771316,9.37157936209166e-07,2.24859318643825e-07,0.00431147499779721,9.19473894481353e-07,2.43835271634213e-07,0.00415749374788092,9.00371909464349e-07,2.62435203588668e-07,0.00400351249796426,8.79881794944362e-07,2.80630535117344e-07,0.00384953124804722,8.58034792546362e-07,2.98393102985964e-07,0.00369554999812978,8.34864932439139e-07,3.15695633292456e-07,0.00354156874821194,8.10407605494330e-07,3.32511325310706e-07,0.00338758749829368,7.84700848466687e-07,3.48814378770544e-07,0.00323360624837500,7.57783983289134e-07,3.64579520681000e-07,0.00307962499845589,7.29698818805091e-07,3.79782583462939e-07,0.00292564374853635,7.00488365588329e-07,3.94399977668995e-07,0.00277166249861638,6.70197946818640e-07,4.08409317404661e-07,0.00261768124869599,6.38873996549661e-07,4.21788842195564e-07,0.00246369999877519,6.06565072879315e-07,4.34518085841281e-07,0.00230971874885398,5.73320747074298e-07,4.46577250994174e-07,0.00215573749893238,5.39192512789093e-07,4.57947917322423e-07,0.00200175624901039,5.04232772895678e-07,4.68612372656151e-07,0.00184777499908804,4.68495639145756e-07,4.78554356093847e-07,0.00169379374916534,4.32036022951904e-07,4.87758349839274e-07,0.00153981249924232,3.94910320563085e-07,4.96210352670001e-07,0.00138583124931899,3.57175613402596e-07,5.03897136830944e-07,0.00123184999939539,3.18890234532733e-07,5.10807047096580e-07,0.00107786874947154,2.80113083479468e-07,5.16929227302355e-07,0.000923887499547456,2.40904070494005e-07,5.22254440074354e-07,0.000769906249623181,2.01323550088309e-07,5.26774267767038e-07,0.000615924999698744,1.61432640354669e-07,5.30481947806707e-07,0.000461943749774174,1.21292778704340e-07,5.33371552961806e-07,0.000307962499849506,8.09659142765045e-08,5.35438837150271e-07,0.000153981249924771,4.05141886189047e-08,5.36680400096056e-07,0,0,5.37094538385944e-07]';
                
                
                %         X_sys_former=X_sys;
                Xmat(:,n)=X_sys;
                [wabc,uabc]=uw_ele(ncon,X_sys);
            end
    end
    
st_analysis_type='linear';
% st_analysis_type='nonlinear';

    
switch st_analysis_type
    case 'linear'
        [k3]=LTBT(sNode_q,Node,ncon,uabc,wabc,Fsys);
%     Xmat(:,n+1)=X_sys;
    case 'nonlinear'
%         wabc=zeros(sncon,3);      uabc=zeros(sncon,3);

          [k3]=NTBT(sNode_q,Node,ncon,uabc,wabc,Fsys);   %nonlinear timoshenko beam theory

        
    results(n+1)=iter;
end
        B=k3\M_sys_con1;
        V_0=Xmat(:,n);
        V_1=B*Xmat(:,n);
        
        omeqa=sqrt((V_1'*V_0)/(V_1'*V_1));
        T_omeqa=2*pi/omeqa;
        k0=3*E*I/L^3;
        m0=ro_t*A*L;
        omeqa0=sqrt(k0/m0);
        T0=2*pi/omeqa0;
        T_T(n)=T_omeqa/T0;
        hjkl=0;
        w_r(n)=wmax/r;
end

% inp_graph(nom_inc_F0)=X_load_incre((sNode_q-1)/2*3+2,n+1);
% inp_graph(nom_inc_F0)=X_load_incre((sNode_q-1)/2*3+2,nom_inc_F0);



% if nom_inc_F0==1
%     X_load_incre(:,nom_inc_F0)=X_sys;
% else
% X_load_incre(:,nom_inc_F0)=X_load_incre(:,nom_inc_F0-1)+X_sys;
% end    %if



hold on
% draw= Xmat((sNode_q-1)/2*3+2,:);
draw= Xmat((sNode_q-1)*3+1,:);

% tvec=0:dt:TT-dt;
% plot(tvec,draw,'DisplayName','axial disp');
% 
% draw1= Xmat((sNode_q-1)*3+2,:);
% 
% % tvec=0:dt:TT;
% plot(tvec,draw1,'DisplayName','lateral disp');

plot(w_r,T_T,'DisplayName','Free vibration');



X_sys((sNode_q-1)/2*3+2)


inp_graph(nom_inc_F0)=X_load_incre((sNode_q-1)/2*3+2,nom_inc_F0);


% end     %for nom_inc_F0
%% Showing the element displacements for the maximum applied load.
NI=NI_F0;
%% Plot
% figure
hold on
plot(inp_graph,1:NI_F0)


for el=1:sncon
    n1=ncon(el,1);
    n2=ncon(el,2);
    %   r1=3*(n1-1)+1:3*(n1-1)+3;
    %   r2=3*(n2-1)+1:3*(n2-1)+3;
    
    r1=3*(n1-1+el-1)+1:3*(n1-1+el-1)+3;   % Note: Just for Unbranched frams
    r2=3*(n2-1+el-1)+1:3*(n2-1+el-1)+3;
    r3=3*(n2-1+el)+1:3*(n2-1+el)+3;
    
    rr=[r1 r2 r3];
    x_ele(el,:)=X_load_incre(rr,NI);
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

a


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

function [Fsys]=F_sys(sNode_q,Node,ncon,f0,q0,F_app)
global A E Ks G I sncon  load_glo_dir

  Fsys=zeros(3*sNode_q,1);

for el=1:sncon
    n1=ncon(el,1);    n2=ncon(el,2);
    xa=Node(n1,1);    xc=Node(n2,1);
    

    F_ele=calculate_F_element(xa,xc,f0,q0,F_app);
    
    % Making K_system and F_system
    r1=3*(n1-1+el-1)+1:3*(n1-1+el-1)+3;   % Note: Just for Unbranched frams
    r2=3*(n2-1+el-1)+1:3*(n2-1+el-1)+3;
    r3=3*(n2-1+el)+1:3*(n2-1+el)+3;
    
    rrr=[r1 r2 r3];
    % T_sys(rrr(:),rrr(:))=T_sys(rrr(:),rrr(:))+T_ele(:,:);
    % R_sys(rrr(:))=R_sys(rrr(:))+R_ele(:);
    
    Fsys(rrr(:))=Fsys(rrr(:))+F_ele(:);
    
    
end

% if xa==0
%     F_ele(1,1)=F0;
% elseif floor(xa/dx)==
Fsys(load_glo_dir)=F_app;






function [Ksys]=K_sys(sNode_q,Node,ncon,uabc,wabc)
global A E Ks G I sncon 

Ksys=zeros(3*sNode_q,3*sNode_q);  
for el=1:sncon
    n1=ncon(el,1);    n2=ncon(el,2);
    xa=Node(n1,1);    xc=Node(n2,1);
    
    wabc_e=wabc(el,:);     uabc_e=uabc(el,:);
    
    % K_ele=calculate_K_element(xa,xb,A,E,Ks,G,I,wab_e);
    % [K_ele T_ele]=calculate_K_T_element(xa,xc,A,E,Ks,G,I,wabc_e,uabc_e);
    [K_ele]=calculate_K__element(xa,xc,wabc_e,uabc_e);


    
    % Making K_system and F_system
    r1=3*(n1-1+el-1)+1:3*(n1-1+el-1)+3;   % Note: Just for Unbranched frams
    r2=3*(n2-1+el-1)+1:3*(n2-1+el-1)+3;
    r3=3*(n2-1+el)+1:3*(n2-1+el)+3;
    
    rrr=[r1 r2 r3];
    % T_sys(rrr(:),rrr(:))=T_sys(rrr(:),rrr(:))+T_ele(:,:);
    % R_sys(rrr(:))=R_sys(rrr(:))+R_ele(:);
    
    Ksys(rrr(:),rrr(:))=Ksys(rrr(:),rrr(:))+K_ele(:,:);

    
    
end




function [M_sys]=Mass_sys(sNode_q,Node,ncon)
global A  I sncon ro_t

M_sys=zeros(3*sNode_q,3*sNode_q);
% F_sys=zeros(3*sNode_q,1);

for el=1:sncon
    n1=ncon(el,1);    n2=ncon(el,2);
    xa=Node(n1,1);    xc=Node(n2,1);
    
    %     wabc_e=wabc(el,:);     uabc_e=uabc(el,:);
    
    % K_ele=calculate_K_element(xa,xb,A,E,Ks,G,I,wab_e);
    % [K_ele T_ele]=calculate_K_T_element(xa,xc,A,E,Ks,G,I,wabc_e,uabc_e);
    %     [K_ele]=calculate_K__element(xa,xc,A,E,Ks,G,I,wabc_e,uabc_e);
    
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
    
    M_sys(rrr(:),rrr(:))=M_sys(rrr(:),rrr(:))+M_ele(:,:);
    
    
end





function [K_ele]=calculate_K__element(xa,xc,wabc_e,uabc_e)
global Ay A E Ks G I
Axx=E*A;
Dxx=E*I;
Sxx=Ks*G*Ay;
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
function [K_ele, T_ele]=calculate_K_T_element(xa,xc,A,E,Ks,G,I,wabc_e,uabc_e)
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
function K_ele=calculate_K_element(xa,xb,A,E,Ks,G,I,wab_e)
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

%% function  calculate_F_element
function F_ele=calculate_F_element(xa,xc,f0,q0,F0)
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
F3=zeros(3,1);    %F3=q0*F_F3;
F_ele([1;4;7])=F1;
F_ele([2;5;8])=F2;
F_ele([3;6;9])=F3;
F_ele=F_ele';     %F_ele=double(F_ele');
% F_ele=vpa(F_ele');

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

%% function input1
function [Node ncon Constraints]=input1()
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
global nsec initial_condition st_analysis_type ...
    Load t G nt dt dx A Ay Az I Pcr nx  Ks TT E ro_t nu  coef_damp max_iter L

initial_condition='given';
initial_condition='not_given';

st_analysis_type='linear';
st_analysis_type='nonlinear';

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
% nsec=40;
nsec=20;
L=6;

Node=zeros(nsec+1,3);
Node(:,1)=0:6/nsec:L;
ncon=zeros(nsec,7);
ncon(:,1)=1:nsec;
ncon(:,2)=2:nsec+1;
c1=1e27;
c2=1e-27;
% Constraints=[1	c2	c1	c1
%     nsec+1	c1	c1	c1];
Constraints=[1	c2	c1	c2
    nsec+1	c1	c1	c2];


Ks=5/6;
Ks=1; % Ks is applied in ay and az.

% TT=1;
TT=0.3;
% TT=0.5;
E= 29e9  ;
ro_t= 2500     ;
nu=0.2;
bbbb=0.5;
hhhh=0.3;
ay=1.22;
az=1.2;

coef_damp=0;
omega=80.88    ;
max_iter=1000;
% Load=-0.5*Pcr*cos(2*omega*t);



A=bbbb*hhhh;

Ay=(1/ay)*A;
Az=(1/az)*A;
A=Ay;
Iz=bbbb*hhhh^3/12;
I=Iz;
Pcr=pi*pi*E*Iz/(L*L)     ;

ct=sqrt(E/ro_t);
% NSection=20;
nx=nsec+1;
dx=L/nsec;
dt=5*dx/ct;
% dt=dx/10/ct;

nt=floor(TT/dt)+1;
G=E/(2*(1+nu));
t=0:dt:TT;

Load=0.5*Pcr*cos(2*omega*t);   % load at the boundary
% Load=0.5*Pcr*cos(omega*t);   % load at the boundary






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
% 
% function [Node, ncon, Constraints]=input7()
% global nsec initial_condition st_analysis_type ...
%     Load t G nt dt dx A Ay Az I Pcr nx  Ks TT E ro_t nu  coef_damp max_iter L load_glo_dir
% 
% % initial_condition='given';
% initial_condition='not_given';
% 
% % st_analysis_type='linear';
% st_analysis_type='nonlinear';
% 
% %[node_number rigidity_u rigidity_w rigidity_teta]
% % Constraints=[1	c2	c1	c2
% % 5	c1	c2	c1];
% 
% force=[];
% 
% %%Node  &   ncon
% poisson_ratio=0 ; %No Pois. Coupl. poisson_ratio=0    Poi. Coup. included    =0.3
% 
% %%ncon
% %3: inner diameter 4:epsilon=raphness   5:E  6:G   7:Rho   8:thickness  9:poisson ratio
% 
% nsec=20;
% L=6;
% 
% Node=zeros(nsec+1,3);
% Node(:,1)=0:6/nsec:L;
% ncon=zeros(nsec,7);
% ncon(:,1)=1:nsec;
% ncon(:,2)=2:nsec+1;
% c1=1e27;
% c2=1e-27;
% % Constraints=[1	c2	c1	c1
% %     nsec+1	c1	c1	c1];
% % Constraints=[1	c2	c1	c2
% %     nsec+1	c1	c1	c2];
% % Constraints=[1	c1	c1	c2
% %     nsec+1	c2	c1	c2];
% Constraints=[1	c1	c1	c1
%     nsec+1	c2	c2	c2];
% Ks=5/6;
% Ks=1; % Ks is applied in ay and az.
% 
% % TT=1;
% % TT=0.05;
% TT=0.5;
% E= 29e9  ;
% ro_t= 2500     ;
% nu=0.2;
% bbbb=0.5;
% hhhh=0.3;
% ay=1.22;
% az=1.2;
% 
% coef_damp=0;
% omega=80.88    ;
% max_iter=1000;
% % Load=-0.5*Pcr*cos(2*omega*t);
% 
% 
% 
% A=bbbb*hhhh;
% 
% Ay=(1/ay)*A;
% Az=(1/az)*A;
% % A=Ay;
% Iz=bbbb*hhhh^3/12;
% I=Iz;
% Pcr=pi*pi*E*Iz/(L*L)     ;
% 
% ct=sqrt(E/ro_t);
% % NSection=20;
% nx=nsec+1;
% dx=L/nsec;
% dt=5*dx/ct;
% % dt=dx/10/ct;
% 
% nt=floor(TT/dt)+1;
% G=E/(2*(1+nu));
% t=0:dt:TT;
% 
% % only works for left or right side loadig
% load_node_dir=[nsec+1 1];
% % load_node_dir=[1 1];   %leads to some problems for intermediate concentrated loads; must be revised.
% 
% 
% % load_glo_dir=1;
% % load_glo_dir=121;
% % load_glo_dir=((load_node_dir(1)-1)*2-2)*3+load_node_dir(2);   %not sure, too hard to check now.
% load_glo_dir=2*3*(-1+load_node_dir(1))+load_node_dir(2); %not sure, too hard to check now.
% 
% 
% Load=-1e7*ones(1,length(t));   % load at the boundary
% 
% % Load=-0.5*Pcr*cos(2*omega*t);   % load at the boundary
% % Load=0.5*Pcr*cos(omega*t);   % load at the boundary
% 
% 
% 
% 
% 
% 
% return



function [Node, ncon, Constraints]=input8()
global nsec initial_condition st_analysis_type ...
    Load t G nt dt dx A Ay Az I Pcr nx  Ks TT E ro_t nu  coef_damp max_iter L load_glo_dir r

% initial_condition='given';
initial_condition='not_given';

st_analysis_type='linear';
% st_analysis_type='nonlinear';

%[node_number rigidity_u rigidity_w rigidity_teta]
% Constraints=[1	c2	c1	c2
% 5	c1	c2	c1];

force=[];

%%Node  &   ncon
poisson_ratio=0 ; %No Pois. Coupl. poisson_ratio=0    Poi. Coup. included    =0.3

%%ncon
%3: inner diameter 4:epsilon=raphness   5:E  6:G   7:Rho   8:thickness  9:poisson ratio


Ks=5/6;
Ks=1; % Ks is applied in ay and az.

% TT=1;
% TT=0.05;
TT=0.5;
E= 70e9  ;
ro_t= 2700     ;
nu=0.2963;

bbbb=6e-3;
hhhh=25e-3;
ay=1.74;
az=1.2;


coef_damp=0;
omega=80.88    ;
max_iter=1000;
% Load=-0.5*Pcr*cos(2*omega*t);



A=bbbb*hhhh;

Ay=(1/ay)*A;
Az=(1/az)*A;
% A=Ay;
Iz=bbbb*hhhh^3/12;
I=Iz;
r=sqrt(I/A);

nsec=20;
L=r*10;
% L=r*20;
% L=r*50;
% L=r*100;

Node=zeros(nsec+1,3);
Node(:,1)=0:L/nsec:L;
ncon=zeros(nsec,7);
ncon(:,1)=1:nsec;
ncon(:,2)=2:nsec+1;

c1=1e27;
c2=1e-27;
% Constraints=[1	c2	c1	c1
%     nsec+1	c1	c1	c1];
% Constraints=[1	c2	c1	c2
%     nsec+1	c1	c1	c2];
Constraints=[1	c1	c1	c2;    
    nsec+1	c2	c1	c2];
% Constraints=[1	c1	c1	c1
%     nsec+1	c2	c2	c2];

Pcr=pi*pi*E*Iz/(L*L)     ;

ct=sqrt(E/ro_t);
% NSection=20;
nx=nsec+1;
dx=L/nsec;
dt=dx/ct;
% dt=5*dx/ct;
% dt=dx/10/ct;
TT=1000*dt;
nt=floor(TT/dt)+1;
G=E/(2*(1+nu));
t=0:dt:TT;

% only works for left or right side loadig
load_node_dir=[nsec+1 1];
% load_node_dir=[1 1];   %leads to some problems for intermediate concentrated loads; must be revised.


% load_glo_dir=1;
% load_glo_dir=121;
% load_glo_dir=((load_node_dir(1)-1)*2-2)*3+load_node_dir(2);   %not sure, too hard to check now.
load_glo_dir=2*3*(-1+load_node_dir(1))+load_node_dir(2); %not sure, too hard to check now.


Load=-1e7*ones(1,length(t));   % load at the boundary

% Load=-0.5*Pcr*cos(2*omega*t);   % load at the boundary
% Load=0.5*Pcr*cos(omega*t);   % load at the boundary
return







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
F_M22=F_M33;
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

function [X,Xd,Xdd] =Newmark(M,C,K,F,dt,Xi,Xdi,Xddi)

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
c6 = dt*(1 - gamma_ ) ;
c7 = dt* gamma_  ;
M=M*1;
Keff = c0*M + c1*C + K ;
Feff = F+ M*(c0*Xi+c2*Xdi+c3*Xddi) ...
    +C*(c1*Xi+c4*Xdi+c5*Xddi) ;

% warning off all

% X=tridiagonal_matrix_algorithm(Keff,Feff);

X=Keff\Feff;


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

function [k3]=LTBT(sNode_q,Node,ncon,uabc,wabc,Fsys)
global Xmat Xdmat Xddmat toler iter max_iter np_deg coef_damp n ...
    M_sys_con1 dt
%     while toler>0.001 && iter<max_iter
%         iter=iter+1;
        %% Assembling
        % T_sys=zeros(3*sNode_q,3*sNode_q);    R_sys=zeros(3*sNode_q,1);
        
        
        
        [Ksys]=K_sys(sNode_q,Node,ncon,uabc,wabc);
        
        [K_sys_con1,F_sys_con1]=K_F_sys_con1(Ksys,Fsys,np_deg,sNode_q);  % apply constraints Boundary conds)
        C=coef_damp*K_sys_con1;
%         X_sys_former=X_sys;
         X_sys_former=Xmat(:,n+1);

        
        f=F_sys_con1;
        k3=K_sys_con1;
        
        %         X_sys=K_sys_con1\F_sys_con1;
        Xi=Xmat(:,n);
        Xdi=Xdmat(:,n);
        Xddi=Xddmat(:,n);
        [X,Xd,Xdd] =Newmark(M_sys_con1,C,k3,f,dt,Xi,Xdi,Xddi);
        Xmat(:,n+1)=X;
        Xdmat(:,n+1)=Xd;
        Xddmat(:,n+1)=Xdd;
        %      Xi=X;
        %      Xdi=Xd;
        %      Xddi=Xdd;
        X_sys=X;
        toler=max(abs((X_sys_former-X_sys)./X_sys));
        
        %         toler=max(abs(X_sys_former-X_sys));
        
        %         Xmat(:,n+1)=X_sys;
        [wabc,uabc]=uw_ele(ncon,X_sys);


function [k3]=NTBT(sNode_q,Node,ncon,uabc,wabc,Fsys)
global Xmat Xdmat Xddmat toler iter max_iter np_deg coef_damp n ...
    M_sys_con1 dt
    while toler>0.001 && iter<max_iter
        iter=iter+1;
        %% Assembling
        % T_sys=zeros(3*sNode_q,3*sNode_q);    R_sys=zeros(3*sNode_q,1);
        
        
        
        [Ksys]=K_sys(sNode_q,Node,ncon,uabc,wabc);
        
        [K_sys_con1,F_sys_con1]=K_F_sys_con1(Ksys,Fsys,np_deg,sNode_q);  % apply constraints Boundary conds)
        C=coef_damp*K_sys_con1;
%         X_sys_former=X_sys;
         X_sys_former=Xmat(:,n+1);

        
        f=F_sys_con1;
        k3=K_sys_con1;
        
        %         X_sys=K_sys_con1\F_sys_con1;
        Xi=Xmat(:,n);
        Xdi=Xdmat(:,n);
        Xddi=Xddmat(:,n);
        [X,Xd,Xdd] =Newmark(M_sys_con1,C,k3,f,dt,Xi,Xdi,Xddi);
        Xmat(:,n+1)=X;
        Xdmat(:,n+1)=Xd;
        Xddmat(:,n+1)=Xdd;
        %      Xi=X;
        %      Xdi=Xd;
        %      Xddi=Xdd;
        X_sys=X;
        toler=max(abs((X_sys_former-X_sys)./X_sys));
        
        %         toler=max(abs(X_sys_former-X_sys));
        
        %         Xmat(:,n+1)=X_sys;
        [wabc,uabc]=uw_ele(ncon,X_sys);
        %         Xdmat(:,n+1)=Xdmat(:,n)+dXd;
        %         Xddmat(:,n+1)=Xddmat(:,n)+dXdd;
        %% calculation of wa, wb which are the only nonlinear parameters
        
        % toler_w=max(abs(wab_former-wab(:,2)));
        % wab_former=wab(:,2);
        if iter==max_iter
            break
        end
    end     %while



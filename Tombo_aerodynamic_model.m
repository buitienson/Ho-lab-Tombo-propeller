function [Td_last, Tr, Md_last, ep, sig, dp, Dgr, Td2, G_w, Dgd, Ex, Gxz, Ix, Ez, Iz] = Tombo_aerodynamic_model(nr, PoPP_2_Xcp, Camberline, Materials, Df, ND)
%% Load materials
% Material Parameters
pois = 0.4;
if Materials == 10 %Dragonskin 10
            strain = [1	2	3	4	5	6	7	8	9	10	12	14	16	18	20];
            E_1_10 = [0.113333333	0.113333333	0.111111111	0.11	0.109333333...
                0.108888889	0.106666667	0.105833333	0.105185185	0.104	0.103333333	0.102857143	0.100416667	0.099259259	0.098333333];
            p = polyfit(strain,E_1_10,3);
            Em = polyval(p,0)*10^6; %Young's mudulus of Dragonskin 10 [Pa]
            Gm = Em/2/(1+pois); %Shear's mudulus of Dragonskin 10 [Pa]   
            
        elseif Materials == 20 %Dragonskin 20
            strain = [1	2	3	4	5	6	7	8	9	10	12	14	16	18	20];
            E_1_20 = [0.426666667	0.426666667	0.426666667	0.425	0.424	0.421111111...
                0.417142857	0.415833333	0.411111111	0.407333333	0.401111111	0.395238095	0.39	0.386296296	0.381333333];
            p = polyfit(strain,E_1_20,3);
            Em = polyval(p,0)*10^6; %Young's mudulus of Dragonskin 20 [Pa]
            Gm = Em/2/(1+pois); %Shear's mudulus of Dragonskin 20 [Pa]
            
        elseif Materials == 30 %DragonSkin 30
            strain = [1	2	3	4	5	6	7	8	9	10	12	14	16	18	20];
            E_1_30 = [0.58	0.58	0.573333333	0.568333333	0.562666667	0.556666667	0.551428571...
                0.545833333	0.540740741	0.54	0.531666667	0.523333333	0.514583333	0.506666667	0.499];
            p = polyfit(strain,E_1_30,3);
            Em = polyval(p,0)*10^6; %Young's mudulus of Dragonskin 30 [Pa]
            Gm = Em/2/(1+pois); %Shear's mudulus of Dragonskin 30 [Pa]
end

%% Airfoil information input
% Environment and operating conditions
phi                 = 0;           %angle of mean flow [degree]
rho                 = 1.2041;      %air density [kg/m^3]
% --- Blade parameters
r                   = 10*10^-3;    %position of start point[m]
R                   = 110*10^-3;   %position of end point[m]
%--------------------------------------------------------------------------------------------------
% --- Nodus parameters
xRN                 = 29*10^-3;    %position of representative section [m]
xN                  = 23*10^-3;    %position of start point of Nodus [m]
XN                  = 35*10^-3;    %position of end point of Nodus [m]
LN                  = XN-xN;       %length of Nodus [m]

% --- Nodust Modeling
%Representative section parameters
AgO                 = (270-124.063)*pi/180; %angle of orginal moment axis [rad]
A_RN                = 65.537*10^-6;         %area of representative sectiton [m^2]
I_x                 = 107.979*10^-12;       %Principal moments of inertia of the area, at the centroid (m^4)
I_y                 = 1376.730*10^-12;      %Principal moments of inertia of the area, at the centroid (m^4)

% Centroid of representative section
x_cos = 29*10^(-3); %mm
y_cos = 3.94*10^(-3); %mm
z_cos = 6.99*10^(-3); %mm

% Center of mass of hard wing
x_com = 66.64*10^(-3); %mm
y_com = 0.79*10^(-3); %mm
z_com = 6.23*10^(-3); %mm

%Tranformation of moment of inertia
tr                  = (270-124.063)*pi/180;                                  %Angle of Transformation rotation [rad]
I_xy                = 0;                                                     %In case of revolution bodies
I_uv                = (I_x-I_y)*sin(2*tr)/2+I_xy*cos(2*tr);                  %The product of inertia
I_u                 = (I_x+I_y)/2+(I_x-I_y)*cos(2*tr)/2-I_xy*sin(2*tr);      %The moment of inertia
I_v                 = (I_x+I_y)/2-(I_x-I_y)*cos(2*tr)/2+I_xy*sin(2*tr);      %The moment of inertia
Ix                  = I_u;                       %[m^4]
Iz                  = I_v;                       %[m^4]    

%Composite propotions
% Df                  = 0.94*10^-3;                %diameter of f_fiber [mm]
Df                  = Df*10^(-3);
Vf                  = (Df^2*pi*ND/4)/A_RN; %Volume fraction of fiber [100%]
Vm                  = 1-Vf;                          %Volume fraction of Matrix materials [100%]

%Fiber materials' parameters: Nilon fiber
E11f                = 2.95*10^9;                 %Young's mudulus of Fiber [Pa]
E22f                = 2.95*10^9;                 %Young's mudulus of Fiber [Pa]
G23f                = E22f/2/(1+0.35);                 %Shear mudulus of Fiber [Pa]
E11                 = Vf*E11f+Vm*Em;                 %
E22                 = Em/(1-(Vf^0.5)*(1-(Em/E22f)));  %Chamis model[Pa]
G23                 = Gm/(1-(Vf^0.5)*(1-Gm/G23f));    %Chamis model [Pa] 
Ex = E22;
Ez = E22;
Gxz = G23;

F_ctl   = 4*nr^2*pi^2*x_com*5.715*10^(-3)/60^2/4; %inertial force
G_w = 5.715*10^(-3)*9.8; %N - Gravity force

%--------------------------------------------------------------------------------------------------
% --- Parameter conversion
omega               = nr*pi*2/60;                 %angular velocity of rotor [rad/s] 
phi_rad             = phi*pi/180;                %[rad]
%--------------------------------------------------------------------------------------------------
%% Thrust force calculation
%% Rigid propeller
%Part 1: Point 1-4, order 2
x1                  = PoPP_2_Xcp(1:4,2);             %Positions of section; 
y11                 = PoPP_2_Xcp(1:4,3);             %Points of leading edge;
y12                 = PoPP_2_Xcp(1:4,5) ;            %Points of trailing edge;
theta1              = PoPP_2_Xcp(1:4,7) + phi_rad;   %Angle of attach;
xcp1                = PoPP_2_Xcp(1:4,9);             %Positions of aerodynamic center;
zcp1                = PoPP_2_Xcp(1:4,8);             %Positions of aerodynamic center;


y11_p               = polyfit(x1,y11,2);         %Call polyfit to generate a cubic fit to predict y11 from x1
y11_fit             = polyval(y11_p,x1);        %Compute the residual values as a vector of signed numbers
y11_resid           = y11 - y11_fit;            %Square the residuals and total them to obtain the residual sum of squares
y11_SSresid         = sum(y11_resid.^2);         %Square the residuals and total them to obtain the residual sum of squares
y11_SStotal         = (length(y11)-1) * var(y11); %Compute the total sum of squares of y1 by multiplying the variance of y11
y11_rsq             = 1 - y11_SSresid/y11_SStotal; %Compute simple R2 for the cubic fit
y11_rsq_adj         = 1 - y11_SSresid/y11_SStotal * (length(y11)-1)/(length(y11)-length(y11_p)); %Compute adjusted R2

y12_p               = polyfit(x1,y12,2);         %Call polyfit to generate a cubic fit to predict y12 from x1
y12_fit             = polyval(y12_p,x1);        %Compute the residual values as a vector of signed numbers
y12_resid           = y12 - y12_fit;            %Square the residuals and total them to obtain the residual sum of squares
y12_SSresid         = sum(y12_resid.^2);         %Square the residuals and total them to obtain the residual sum of squares
y12_SStotal         = (length(y12)-1) * var(y12); %Compute the total sum of squares of y1 by multiplying the variance of y12
y12_rsq             = 1 - y12_SSresid/y12_SStotal; %Compute simple R2 for the cubic fit
y12_rsq_adj         = 1 - y12_SSresid/y12_SStotal * (length(y12)-1)/(length(y12)-length(y12_p)); %Compute adjusted R2

theta1_p            = polyfit(x1,theta1,2);         %Call polyfit to generate a cubic fit to predict theta1 from x1
theta1_fit          = polyval(theta1_p,x1);        %Compute the residual values as a vector of signed numbers
theta1_resid        = theta1 - theta1_fit;         %Square the residuals and total them to obtain the residual sum of squares
theta1_SSresid      = sum(theta1_resid.^2);         %Square the residuals and total them to obtain the residual sum of squares
theta1_SStotal      = (length(theta1)-1) * var(theta1); %Compute the total sum of squares of y1 by multiplying the variance of theta1
theta1_rsq          = 1 - theta1_SSresid/theta1_SStotal;  %Compute simple R2 for the cubic fit
theta1_rsq_adj      = 1 - theta1_SSresid/theta1_SStotal * (length(theta1)-1)/(length(theta1)-length(theta1_p)); %Compute adjusted R2

xcp1_p              = polyfit(x1,xcp1,2);         %Call polyfit to generate a cubic fit to predict xcp1 from x1
xcp1_fit            = polyval(xcp1_p,x1);        %Compute the residual values as a vector of signed numbers
xcp1_resid          = xcp1 - xcp1_fit;            %Square the residuals and total them to obtain the residual sum of squares
xcp1_SSresid        = sum(xcp1_resid.^2);         %Square the residuals and total them to obtain the residual sum of squares
xcp1_SStotal        = (length(xcp1)-1) * var(xcp1);%Compute the total sum of squares of y1 by multiplying the variance of xcp1
xcp1_rsq            = 1 - xcp1_SSresid/xcp1_SStotal; %Compute simple R2 for the cubic fit
xcp1_rsq_adj        = 1 - xcp1_SSresid/xcp1_SStotal * (length(xcp1)-1)/(length(xcp1)-length(xcp1_p)); %Compute adjusted R2

zcp1_p              = polyfit(x1,zcp1,2);         %Call polyfit to generate a cubic fit to predict zcp1 from x1
zcp1_fit            = polyval(zcp1_p,x1);        %Compute the residual values as a vector of signed numbers
zcp1_resid          = zcp1 - zcp1_fit;            %Square the residuals and total them to obtain the residual sum of squares
zcp1_SSresid        = sum(zcp1_resid.^2);         %Square the residuals and total them to obtain the residual sum of squares
zcp1_SStotal        = (length(zcp1)-1) * var(zcp1);%Compute the total sum of squares of y1 by multiplying the variance of zcp1
zcp1_rsq            = 1 - zcp1_SSresid/zcp1_SStotal; %Compute simple R2 for the cubic fit
zcp1_rsq_adj        = 1 - zcp1_SSresid/zcp1_SStotal * (length(zcp1)-1)/(length(zcp1)-length(zcp1_p));%Compute adjusted R2

gfbeta1_0           = y11_p(1)- y12_p(1);
gfbeta1_1           = y11_p(2)- y12_p(2);
gfbeta1_2           = y11_p(3)- y12_p(3);

T1                  = @(t) 1.2041*omega^2*(gfbeta1_2 + gfbeta1_1.*t + gfbeta1_0.*t.^2).*...
    (0.225 + 1.58.*sin((2.13.*(theta1_p(3) + theta1_p(2).*t + theta1_p(1).*t.^2)-7.2).*(pi/180))).*(t.^2);
Dg1                 = @(t) 1.2041*omega^2*(gfbeta1_2 + gfbeta1_1.*t + gfbeta1_0.*t.^2).*...
    (1.92 - 1.55.*cos((2.04.*(theta1_p(3) + theta1_p(2).*t + theta1_p(1).*t.^2)-9.82).*(pi/180))).*(t.^2);
%Part 2: Point 4-24, order 3
x2                  = PoPP_2_Xcp(4:24,2);                %Positions of section [m]; 
y21                 = PoPP_2_Xcp(4:24,3);                %Points of leading edge [m];
y22                 = PoPP_2_Xcp(4:24,5);                %Points of trailing edge [m];
theta2              = PoPP_2_Xcp(4:24,7);                %twist angle of airfoil [degree];
xcp2                = PoPP_2_Xcp(4:24,8);                %Positions of aerodynamic center;
zcp2                = PoPP_2_Xcp(4:24,9);                %Positions of aerodynamic center;
chord2              = PoPP_2_Xcp(4:24,10);               %chord length;

y21_p               = polyfit(x2,y21,3) ;                %Call polyfit to generate a cubic fit to predict y21 from x2
y21_fit             = polyval(y21_p,x2);                %Compute the residual values as a vector of signed numbers
y21_resid           = y21 - y21_fit;                    %Square the residuals and total them to obtain the residual sum of squares
y21_SSresid         = sum(y21_resid.^2)    ;             %Square the residuals and total them to obtain the residual sum of squares
y21_SStotal         = (length(y21)-1) * var(y21) ;       %Compute the total sum of squares of y1 by multiplying the variance of y21
y21_rsq             = 1 - y21_SSresid/y21_SStotal ;      %Compute simple R2 for the cubic fit
y21_rsq_adj         = 1 - y21_SSresid/y21_SStotal * (length(y21)-1)/(length(y21)-length(y21_p));%Compute adjusted R2

y22_p               = polyfit(x2,y22,3);                %Call polyfit to generate a cubic fit to predict y22 from x2
y22_fit             = polyval(y22_p,x2);                %Compute the residual values as a vector of signed numbers
y22_resid           = y22 - y22_fit;                    %Square the residuals and total them to obtain the residual sum of squares
y22_SSresid         = sum(y22_resid.^2)  ;               %Square the residuals and total them to obtain the residual sum of squares
y22_SStotal         = (length(y22)-1) * var(y22) ;       %Compute the total sum of squares of y1 by multiplying the variance of y22
y22_rsq             = 1 - y22_SSresid/y22_SStotal;       %Compute simple R2 for the cubic fit
y22_rsq_adj         = 1 - y22_SSresid/y22_SStotal * (length(y22)-1)/(length(y22)-length(y22_p));%Compute adjusted R2

theta2_p            = polyfit(x2,theta2,3);              %Call polyfit to generate a cubic fit to predict theta2 from x2
theta2_fit          = polyval(theta2_p,x2);             %Compute the residual values as a vector of signed numbers
theta2_resid        = theta2 - theta2_fit;              %Square the residuals and total them to obtain the residual sum of squares
theta2_SSresid      = sum(theta2_resid.^2);              %Square the residuals and total them to obtain the residual sum of squares
theta2_SStotal      = (length(theta2)-1) * var(theta2);  %Compute the total sum of squares of y1 by multiplying the variance of theta2
theta2_rsq          = 1 - theta2_SSresid/theta2_SStotal ;%Compute simple R2 for the cubic fit
theta2_rsq_adj      = 1 - theta2_SSresid/theta2_SStotal * (length(theta2)-1)/(length(theta2)-length(theta2_p)); %Compute adjusted R2

xcp2_p              = polyfit(x2,xcp2,3);                %Call polyfit to generate a cubic fit to predict xcp2 from x2
xcp2_fit            = polyval(xcp2_p,x2);               %Compute the residual values as a vector of signed numbers
xcp2_resid          = xcp2 - xcp2_fit;                  %Square the residuals and total them to obtain the residual sum of squares
xcp2_SSresid        = sum(xcp2_resid.^2);                %Square the residuals and total them to obtain the residual sum of squares
xcp2_SStotal        = (length(xcp2)-1) * var(xcp2);      %Compute the total sum of squares of y1 by multiplying the variance of xcp2
xcp2_rsq            = 1 - xcp2_SSresid/xcp2_SStotal;     %Compute simple R2 for the cubic fit
xcp2_rsq_adj        = 1 - xcp2_SSresid/xcp2_SStotal * (length(xcp2)-1)/(length(xcp2)-length(xcp2_p));%Compute adjusted R2

zcp2_p              = polyfit(x2,zcp2,3);               %Call polyfit to generate a cubic fit to predict zcp2 from x2
zcp2_fit            = polyval(zcp2_p,x2);               %Compute the residual values as a vector of signed numbers
zcp2_resid          = zcp2 - zcp2_fit;                  %Square the residuals and total them to obtain the residual sum of squares
zcp2_SSresid        = sum(zcp2_resid.^2) ;               %Square the residuals and total them to obtain the residual sum of squares
zcp2_SStotal        = (length(zcp2)-1) * var(zcp2) ;     %Compute the total sum of squares of y1 by multiplying the variance of zcp2
zcp2_rsq            = 1 - zcp2_SSresid/zcp2_SStotal ;   %Compute simple R2 for the cubic fit
zcp2_rsq_adj        = 1 - zcp2_SSresid/zcp2_SStotal * (length(zcp2)-1)/(length(zcp2)-length(zcp2_p));%Compute adjusted R2

chord2_p            = polyfit(x2,chord2,3)  ;            %Call polyfit to generate a cubic fit to predict chord2 from x2
chord2_fit          = polyval(chord2_p,x2);             %Compute the residual values as a vector of signed numbers
chord2_resid        = chord2 - chord2_fit;              %Square the residuals and total them to obtain the residual sum of squares
chord2_SSresid      = sum(chord2_resid.^2) ;             %Square the residuals and total them to obtain the residual sum of squares
chord2_SStotal      = (length(chord2)-1) * var(chord2) ; %Compute the total sum of squares of y1 by multiplying the variance of chord2
chord2_rsq          = 1 - chord2_SSresid/chord2_SStotal; %Compute simple R2 for the cubic fit
chord2_rsq_adj      = 1 - chord2_SSresid/chord2_SStotal * (length(chord2)-1)/(length(chord2)-length(chord2_p)); %Compute adjusted R2

%% Define Center of pressure
X_cop = PoPP_2_Xcp(4:24,11); 
Z_cop = PoPP_2_Xcp(4:24,12);
Y_cop = PoPP_2_Xcp(4:24,13);
X_cop_p            = polyfit(x2,X_cop,3);
Z_cop_p            = polyfit(x2,Z_cop,3);
Y_cop_p            = polyfit(x2,Y_cop,3);

%% Define f(x)-g(x)
gfbeta2_0           = y21_p(1)- y22_p(1);
gfbeta2_1           = y21_p(2)- y22_p(2);
gfbeta2_2           = y21_p(3)- y22_p(3);
gfbeta2_3           = y21_p(4)- y22_p(4);

T2                  = @(t) 1.2041*omega^2*(gfbeta2_3 + gfbeta2_2.*t + gfbeta2_1.*t.^2+ gfbeta2_0.*t.^3).*...
    (0.225 + 1.58.*sin((2.13.*(theta2_p(4) + theta2_p(3).*t + theta2_p(2).*t.^2+ theta2_p(1).*t.^3)-7.2).*...
    (pi/180))).*(t.^2);
Dg2                 = @(t) 1.2041*omega^2*(gfbeta2_3 + gfbeta2_2.*t + gfbeta2_1.*t.^2+ gfbeta2_0.*t.^3).*...
    (1.92 - 1.55.*cos((2.04.*(theta2_p(4) + theta2_p(3).*t + theta2_p(2).*t.^2+ theta2_p(1).*t.^3)-9.82).*...
    (pi/180))).*(t.^2);

%% Total thrust force of rigid propeller
Tr1                 = integral(T1,x1(1),x1(4));
Tr2                 = integral(T2,x2(1),x2(21));
Tr                  = Tr1+Tr2  ;                 %Total thrust force of rigid propeller
Dgr1                = integral(Dg1,x1(1),x1(4));
Dgr2                = integral(Dg2,x2(1),x2(21));
Dgr                 = Dgr1+Dgr2  ;               %Total drag force of rigid propeller

%% Deformable propeller
%Part 1: Point 1-4
Td1                 = Tr1+integral(T2,x2(1),XN);
Dgd1                = Dgr1;

%Part 2: Point 4-24
%% Compute ep and sig
ep = atan(y_com/x_com)*0.8;
sig = atan(z_com-4.5*10^(-3))/x_com*0.8;
ep_g=ep;
sig_g=sig;

%% Parameters of Camberline
% Define F
tc          = Camberline(1:4,2);   %thickness normal to camberline                            
uc          = Camberline(1:4,3);   %elementary length along camber line                              
tc2         = Camberline(4:11,2);  %thickness normal to camberline                            
uc2         = Camberline(4:11,3);  %elementary length along camber line 
A           = 65.57*10^(-6);       %area of section
%Build t(u)
tc_p        = polyfit(uc,tc,2);
tc_fit      = polyval(tc_p,uc);
tc_resid    = tc - tc_fit;
tc_SSresid  = sum(tc_resid.^2);
tc_SStotal  = (length(tc)-1) * var(tc);
tc_rsq      = 1 - tc_SSresid/tc_SStotal;
tc_rsq_adj  = 1 - tc_SSresid/tc_SStotal * (length(tc)-1)/(length(tc)-length(tc_p));
tc2_p       = polyfit(uc2,tc2,3);
F1          = @(u) (tc(3) + tc(2).*u.^1 + tc(1).*u.^2).^3;
F2          = @(u) (tc2_p(4) + tc2_p(3).*u + tc2_p(2).*u.^2 + tc2_p(1).*u.^3).^3;
I_F1        = integral(F1,uc(4),uc(1));
I_F2        = integral(F2,uc2(8),uc2(1));
F           = I_F1 + I_F2;

%Define U
U           = 36.58*10^-3; %length of camberline [m]    

%%     --- Step 1: Initial guess
Tor         = Tr2*5*10^(-3);
dp          = 3*(1+4*F/(3*A*U^2))*Tor*LN/(Gxz*F);
dp_g        = dp;
Td2 = Tr2;
Dgd2 = Dgr2;
CoS = [x_cos y_cos z_cos]';
CoM_design = [x_com y_com z_com]';
CoM = eul2rotm([0 -0.0271 0]) * (CoM_design - CoS) + CoS;

%%     --- Step 2: Loops
%Algorithm options
nbIt                = 1000; % Maximum number of iterations
Tol                 = 0.01;
for i=1:nbIt
    %% Compute deformable angle
    syms t
    
    % --- Step 1: Compute sigma 
    eqn_sig = t - (-F_ctl*tan(t - 0.095 - 0.0271)+ Td2 - G_w)...
        *(XN/x_com)*LN^2/(2*Ex*Ix);
    AA = (XN/x_com)*LN^2/(2*Ex*Ix);
    V_sig = vpasolve(eqn_sig,t,[-1.57 - (- 0.095 - 0.0271), 1.57 - (- 0.095 - 0.0271)]);
    sig = double(V_sig);

    % --- Step 2: Compute epsilon
    eqn_ep = t - (F_ctl * tan(0.1418-t)+ Dgd2)*(x_com/XN)*LN^2/(2*Ez*Iz);
    BB =(x_com/XN)*LN^2/(2*Ez*Iz)  ;
    V_ep = vpasolve(eqn_ep,t,[-1.57 1.57]);
    ep = double(V_ep);
  
    % --- Step 3: Compute dp
    F1 = (-F_ctl*tan(sig - 0.095 - 0.0271)+ Td2-G_w);
    F2 = (F_ctl*tan(ep)+ Dgd2);
    Tor = F1*(y_cos - y_com) + F2*(z_cos - z_com);
    dp = 3*(1+4*F/(3*A*U^2))*Tor*LN/(Gxz*F);
    
    %% update CoM
    CoM = eul2rotm([ep sig dp]) * (CoM - CoS) + CoS;
    
    %% update XN
    Nodus = [XN 0 0]' ;
    Nodus = eul2rotm([ep sig dp]) * (Nodus - CoS) + CoS;
    XN = Nodus(1);
    
    %% Compute Td, Dfd  
    Theta_fun_deg = @(u) (theta2_p(4) + theta2_p(3).*u...
        + theta2_p(2).*u.^2+ theta2_p(1).*u.^3);
    Theta_fun   = @(u) Theta_fun_deg(u)*pi/180; %radian

    % Compute Td
    Td2_f1 = @(u) (gfbeta2_3 + gfbeta2_2.*u + gfbeta2_1.*u.^2+ gfbeta2_0.*u.^3);
    Td2_f2 = @(u) cos(Theta_fun(u)- dp);
    Td2_f3 = @(u) cos(Theta_fun(u));
    CL_fun = @(u) 0.225 + 1.58.*sin((2.13.*Theta_fun_deg(u)-7.2).*pi/180-2.13*dp);        
    Td2_fun    = @(u) rho*omega^2.*Td2_f1(u).* ...
        cos(sig).*cos(ep).*Td2_f2(u)./Td2_f3(u)...
            .*CL_fun(u).*(u.^2);   
    Td2         = integral(Td2_fun,XN,R);
    Td          = Td1+Td2;
    
    %% Compute Dgd
    Dgd2_f1 = @(u) (gfbeta2_3 + gfbeta2_2.*u + gfbeta2_1.*u.^2+ gfbeta2_0.*u.^3);
    Dgd2_f2 = @(u) cos(Theta_fun(u)+ dp);
    Dgd2_f3 = @(u) cos(Theta_fun(u));
    CD_fun = @(u) 1.92 - 1.55.*cos((2.04.*(Theta_fun_deg(u))-9.82).*pi/180-2.04*dp);   
    Dgd2_fun    = @(u) rho*omega^2*...
    cos(sig)*cos(ep).*Dgd2_f1(u).* Dgd2_f2(u)./Dgd2_f3(u)...
    .*CD_fun(u).*(u.^2);  
    Dgd2        = integral(Dgd2_fun,XN,R);       
    Dgd         = Dgd1 + Dgd2;

%%     --- Step 3: Storing last values
    Td1_last = Td1;
    Td_last     = Td;
    Dgd_last    = Dgd;
    TDratio     = Td/Dgd;
    Md_last     = Tor;
    Ep_last     = ep;
    Sig_last    = sig;
    Dp_last     = dp;    
    
    % --- Convergence Criteria
    if (i>3 && abs((Td-Td_last)/Td) + abs((Dgd-Dgd_last)/Dgd)<Tol); break; end
    if(i==nbIt); fprintf('Maximum iterations reached at Td = ',Td_last); end
    
end %iterative loop for one element

end %Total thrust force of deformable propeller
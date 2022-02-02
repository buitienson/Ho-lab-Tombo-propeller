clc;
clf;
%% Load data
% Airfoil data
PoPP_table = readtable('Airfoil.csv', 'NumHeaderLines', 1);
PoPP = PoPP_table{:, :};

% Representative cross-section data
load Camberline.dat;

%% INPUT
% Matrix materials were denoted follow: DragonSkin 10 as 10, DragonSkin 20
% as 20, and DragonSkin 30 as 30. Note that the length and position of
% nodus are fixed.
matrix_material = 10;%

% Reinforced fiber diameter
fiber_diameter = 0.94; % mm

% Number of tendon fibers
fiber_number = 5;

% Rotational speed of propeller
rotational_speed = 4000; % rpm

%% OUTPUT
[Td_last, Tr, Md_last, ep, sig, dp, Dgr, Td2, G_w, Dgd, Ex, Gxz]...
    = Tombo_aerodynamic_model(rotational_speed + 0.01, PoPP,...
    Camberline, matrix_material, fiber_diameter, fiber_number);

C_Ep = ep * 180/pi;
C_Sig = sig * 180/pi;
C_Dp = dp * 180/pi;

L_over = Td_last/Dgd;
disp(strcat('Thrust force of Tombo propeller: ', num2str(Td_last), ' [N]'))
disp(strcat('Thrust force of a rigid propeller: ', num2str(Tr), ' [N]'))
disp(strcat('Torque on the wing: ', num2str(Md_last), ' [Nm]'))
disp(strcat('Deformable angle alpha: ', num2str(C_Ep), ' [degree]'))
disp(strcat('Deformable angle beta: ', num2str(C_Sig), ' [degree]'))
disp(strcat('Deformable angle sigma: ', num2str(C_Dp), ' [degree]'))
disp(strcat('Young modulus of nodus: ', num2str(Ex/1e6), ' [MPa]'))
disp(strcat('Shear modulus of nodus: ', num2str(Gxz/1e6), ' [MPa]'))
disp(strcat('Lift-over-drag ratio: ', num2str(L_over)))

%% E_field Error Comparison %%
clear all;
close all;


%%%%%%%%%% Load mat files from Python code and obtain fields %%%%%%%%%
load('E_farfield_matrix_sphere_z_offset_test10.mat');
Eff_new = E_farfield_matrix_sphere_z_offset(1,1:180);
load('E_farfield_matrix_sphere_z_offset_old.mat');
Eff_old = E_farfield(1,1:180);

%%%%%%%%%% Define theta as in Python code %%%%%%%%%
theta = linspace(0.001,2*pi,180);

Eff_diff = abs(Eff_new-Eff_old);

figure(1)
plot(theta,abs(Eff_new),theta,abs(Eff_old))
xlabel('$\theta$ [rad]','Interpreter','latex')
legend('$|E_{new}|$ [V/m]','$|E_{old}|$ [V/m]','interpreter','latex')
grid on

figure(2)
plot(theta,Eff_diff)
xlabel('$\theta$ [$^{\circ}$]','Interpreter','latex')
grid on
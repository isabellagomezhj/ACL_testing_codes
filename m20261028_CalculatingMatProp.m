%% Calculating Material Properties
% This code calculates the material properties for the PCL test
% Fill in the directory for the file being used. This code works
% specifically for the formatting in : "S:\BiomechanicsResearch\groupImhauser\ACL Structure Mechanics\Data_Reduced\pcl_Distances_Strains.xlsx"
% and will not work for another file that is formatted differently or has
% different sheets

close all
clear 
clc

filename = "S:\BiomechanicsResearch\groupImhauser\ACL Structure Mechanics\Data_Reduced\pcl_Distances_Strains.xlsx";
sheet = sheetnames(filename);

fig_count=1;
%% Geometric Data
%cross sectional area should be calculated using "CalculatingCSA.m" with
%.stl files collected using the 3D Scanner

CSA = 56.8236953460548;

% %% Ramp - Manual
% % Material properties for the ramp are calculated using the equation :
% % E=stress/strain
% 
% load sheets of interest and remove NaN rows
sheet_7=readtable(filename,'Sheet',sheet(7));
sheet_7(1,:)=[];

% % extract data of interest
% f_acumen=sheet_7.AxialForce;
% r_strain=(sheet_7.strainy+sheet_7.strainy_1+sheet_7.strainy_2)./3;
time=rmmissing(sheet_7.attribute);

% % calculate stress in MPa and resample to fit length of DIC data
% r_stress_full=f_acumen./CSA;
% ds_to=1:length(r_strain);
% ds_from=1:length(r_stress_full);
% r_stress=transpose(interp1(ds_from,r_stress_full,ds_to));
% 
% % calculate young's modulus - will not work, Acumen and DIC output
% % different amounts of data
% E_ramp=r_stress./r_strain;
% 
% % calculate line of best fit and plot both
% figure(fig_count)
% plot(time,E_ramp,'x');
% xlabel('Time [s]');
% ylabel('E-Modulus [MPa]');
% title("Young's Modulus - Ramp");
% 
% figure(fig_count+1)
% fit_2=fitlm(r_strain,r_stress);
% plot(fit_2);
% xlabel('Strain [%]');
% ylabel('Stress [MPa]');
% title("Strain vs. Stress");

%% QLV Method
%References : "S:\BiomechanicsResearch\groupImhauser\ACL Structure Mechanics\Literature\Abramowitch FiniteStrainRate StressRelaxation 2004.pdf"

%%read in relevant sheets%%
sheet_5=readtable(filename,'Sheet',sheet(5));
sheet_5(1,:)=[];
sheet_8=readtable(filename,'Sheet',sheet(8));
sheet_16=readtable(filename,'Sheet',sheet(16));

%%establish time ranges%%
t=sheet_16.TimeZeroed;
N=sheet_16.N;

%find the local min from time 0 to time 2
mask1=(t>=0)&(t<=2);
t_1=t(mask1);
N_1=N(mask1);
min1=min(N_1);
idx1=find(N==min1);

%find the local max from time 2 to time 5
mask2=(t>=2)&(t<=5);
t_2=t(mask2);
N_2=N(mask2);
max1=max(N_2);
idx2=find(N==max1);

%define ramp time and forces
t_ramp=t(idx1:idx2);
N_ramp=N(idx1:idx2);

%find the local max from time 32 to time 35
mask3=(t>=32)&(t<=35);
t_3=t(mask3);
N_3=N(mask3);
max2=max(N_3);
idx3=find(N==max2);

%define stress relaxation time and force
t_relax=t(idx2:idx3);
N_relax=N(idx2:idx3);

%define any critical points
t_ra=t_ramp-t_ramp(1);
t_rx=t_relax-t_ramp(1);
t_0=t_rx(1);

%plot (for verification)
figure(fig_count+1)
plot(t,N,t_1(idx1),N_1(idx1),'r*');
hold on
plot(t(idx2),N(idx2),'r*');
plot(t_ramp,N_ramp,'g');
plot(t(idx3),N(idx3),'r*');
plot(t_relax,N_relax,'m');

%%calculate strain rate%%
avg_str=mean(rmmissing(sheet_8.avgStrainy));
gam=avg_str/mean(rmmissing(sheet_8.time_DIC));

%%define initial guesses%%
x0=[.1 .1 .1 .1 .1];

%%define the integral terms%%
intgnd=@(t,tau,tau_1,tau_2,B,C,gamma)...
    (1+C*(expint((t-tau)/tau_2)-expint((t-tau)/tau_1))).*exp(B*gamma*tau);

%%define eqn5 and eqn6%%
%Ramp
eqn5=@(A,B,C,tau_1,tau_2) arrayfun(@(t_step)...
    ((A*B*gam)/(1+C*log(tau_2/tau_1)))*...
    integral(@(tau) intgnd(t_step,tau,tau_1,tau_2,B,C,gam),0,t_step),t_ra);

R=N_ramp/CSA;

%Stress Relaxation
eqn6=@(A,B,C,tau_1,tau_2) arrayfun(@(t_step)...
    ((A*B*gam)/(1+C*log(tau_2/tau_1)))*...
    integral(@(tau) intgnd(t_step,tau,tau_1,tau_2,B,C,gam),0,t_0),t_rx);

S=N_relax/CSA;

%%nonlinear least-squares regression%%
%levenberg-marquardt was the model used by Abramowitch et al.
min_fun=@(x)...
    [R-eqn5(x(1),x(2),x(3),x(4),x(5));...
    S-eqn6(x(1),x(2),x(3),x(4),x(5))];

options=optimoptions('lsqnonlin','Algorithm','levenberg-marquardt',...
    'MaxIterations',4.00e+3);
[vars,resnorm]=lsqnonlin(min_fun,x0,[],[],options);

%plot for verification
figure(fig_count+1)
plot(t_ramp,R,'.b');
hold on
plot(t_relax,S,'.r');

figure(fig_count+1)
plot(t_ramp,eqn5(vars(1),vars(2),vars(3),vars(4),vars(5)));
hold on
plot(t_ramp,R,'.b');

figure(fig_count+1)
plot(t_relax,eqn6(vars(1),vars(2),vars(3),vars(4),vars(5)));
hold on
plot(t_relax,S,'.r');

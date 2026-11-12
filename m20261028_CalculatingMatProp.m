%% Simplified Quasi-Linear Viscoelastic Theory (QLV) 
% This code utilizes the QLV method found in reference [1] to estimate the
% material properties of a ligament 
%
%
% References : 
% [1] "S:\BiomechanicsResearch\groupImhauser\ACL Structure Mechanics\Literature\Abramowitch FiniteStrainRate StressRelaxation 2004.pdf"
%
%
% Fill in the directory for the file being used. This code works
% specifically for files formatted the same as:
% "S:\BiomechanicsResearch\groupImhauser\ACL Structure Mechanics\Data_Reduced\QLVInputs_07012026.xlsx"
% and will not work for another file that is formatted differently

close all
clear 
clc

filename = "S:\BiomechanicsResearch\groupImhauser\ACL Structure Mechanics\Data_Reduced\QLVInputs_07012026.xlsx";
sheet = sheetnames(filename);

%% Geometric Data
% cross sectional area should be calculated using "CalculatingCSA.m" with
% .stl files collected using the 3D Scanner
CSA = 56.8236953460548;

%% Acquire relevant data from file
sheet_1=readtable(filename,'Sheet',sheet(1)); %DIC data (raw)
sheet_2=readtable(filename,'Sheet',sheet(2)); %Acumen data (filtered)

%% Establish relevant time ranges
% The time ranges for the DIC and Acumen are determined differently. The
% Acumen utilizes a manual masking method to detect the points of interest
% while the DIC utilizes MATLAB's ischange function on the smoothed data.

%%Acumen%%
% the Acumen data must be cut to include just before the ramp begins
% (before the dip that starts the ramp) and to just after the stress
% relaxation ends (include the decrease that occurs at the end of the
% relaxation plateau). Ensure the time vector for the cut data begins at
% t=0.

ACU_t=sheet_2.TimeZeroed; %location of filtered and zeroed Acumen time
N=sheet_2.N; %location of filtered Acumen force outputs

ystrain=rmmissing(sheet_1.AverageStrainY);
DIC_t=linspace(0,4126/25,length(ystrain)); %calculates time based on 25Hz sample rate

% find the local min from time 0 to time 2
mask1=(ACU_t>=0)&(ACU_t<=2);
N1=N(mask1);
idx1=find(N==min(N1));

% find the local max from time 2 to time 5
mask2=(ACU_t>=2)&(ACU_t<=5);
N2=N(mask2);
idx2=find(N==max(N2));

% define ramp time and forces
ACU_t_ramp=ACU_t(idx1:idx2);
N_ramp=N(idx1:idx2);

% find the local max from time 32 to time 35
mask3=(ACU_t>=32)&(ACU_t<=35);
N3=N(mask3);
idx3=find(N==max(N3));

% define stress relaxation time and force
ACU_t_relax=ACU_t(idx2:idx3);
N_relax=N(idx2:idx3);

% filter out overshoot
indexlist=[];
j=1;

for i=1:length(N_ramp)
    if ACU_t_ramp(i)<1.932 && N_ramp(i)>15
        indexlist(j)=i;
        j=j+1;
    end
end

ACU_t_ramp(indexlist)=[];
N_ramp(indexlist)=[];

N_ramp=rmmissing(N_ramp);
ACU_t_ramp=rmmissing(ACU_t_ramp);

% plot to verify overshoot
figure()
plot(ACU_t_ramp,N_ramp,'.');
hold on
plot(ACU_t_relax,N_relax,'.');

% plot (for verification)
figure
hold on
plot(ACU_t,N,ACU_t(idx1),N(idx1),'r*');
plot(ACU_t(idx2),N(idx2),'r*');
plot(ACU_t_ramp,N_ramp,'g');
plot(ACU_t(idx3),N(idx3),'r*');
plot(ACU_t_relax,N_relax,'m');
title('Filtered Data (Ramp & Relaxation) - ACU');

% define any critical points for the regression
t_ra=ACU_t_ramp-ACU_t_ramp(1); %zeroes the ramp time 
t_rx=ACU_t_relax-ACU_t_ramp(1); %matches the relax time to new ramp start
t_0=t_rx(1); %defines t_0 as the relax start (matches ramp end)

%%DIC%%
% The code automatically cuts the DIC data for the test on 07/01/2026.
% Future versions may have to better handle time data to automate how the
% data is read in and filtered to include only relevant points.Ensure the 
% time vector for the cut data begins at t=0.

% cut the data to one cycle
mask4=(DIC_t>=100)&(DIC_t<=132);
DIC_t1=DIC_t(mask4);
DIC_t1=DIC_t1-DIC_t1(1); %zeroes the filtered data
ystrain1=ystrain(mask4);

% plot to verify adequate filtering
figure()
plot(DIC_t1,ystrain1);
title('Filtered Data - DIC');

% smooth data and detect changes
ystrain1_smooth=smoothdata(ystrain1,'movmean');
TF=ischange(ystrain1_smooth,'linear');
idxs=find(TF);

% plot to verify smoothing and critical point identification
figure()
hold on
plot(DIC_t1,ystrain1_smooth);
plot(DIC_t1(TF),ystrain1(TF),'r*');
title('Data Smoothed (W/ Critical Pts.) - DIC');

% define ramp regions (used to calculate gamma)
DIC_t_ramp=DIC_t1(idxs(1):idxs(2));
ystrain_ramp=ystrain1(idxs(1):idxs(2)).*0.01;

% define relax regions
DIC_t_relax=DIC_t1(idxs(2):idxs(3));
ystrain_relax=ystrain1(idxs(2):idxs(3)).*0.01;

% plot to verify sectioning
figure()
hold on
plot(DIC_t_ramp,ystrain_ramp,'m');
plot(DIC_t_relax,ystrain_relax,'g');
plot(DIC_t1(TF),ystrain1(TF)*0.01,'y*');
title('Filtered Data (Ramp & Relaxation) - DIC');

% filter out overshoot
indexlist2=[];
j2=1;

for i=1:length(ystrain_ramp)
    if DIC_t_ramp(i)<5 && ystrain_ramp(i)>0.04
        indexlist2(j2)=i;
        j2=j2+1;
    end
end

DIC_t_ramp(indexlist2)=[];
ystrain_ramp(indexlist2)=[];

ystrain_ramp=rmmissing(ystrain_ramp);
DIC_t_ramp=rmmissing(DIC_t_ramp);

% plot to verify exclusion of overshoot
figure()
plot(DIC_t_ramp,ystrain_ramp);
hold on
plot(DIC_t_relax,ystrain_relax);

%% Calculate strain rate (gamma)
coeffs=polyfit(DIC_t_ramp,ystrain_ramp,1); %fit ramp data (from DIC) linearly
gamma_fit=polyval(coeffs,DIC_t_ramp);
gamma=coeffs(1); %slope of best fit line

% % plot to view fit
% figure()
% hold on
% plot(DIC_t_ra,ystrain_ra,'.');
% plot(DIC_t_ra,gamma_fit);
% title('Gamma');

%% Non-Linear Regression (Levenberg-Marquardt)

% define initial guesses
x0=[0.323937474	3.239374744	32.39374744	0.647874949	64.78749489];

% define the delta_E1 and integral terms
% delta_E1 avoids integral() rounding to a zero numerator which results in
% Inf or NaN values in the regression by evaluating small values at the
% finite limit of the difference (E1()-E1()) found in the integrand
delta_E1=@(x,tau_1,tau_2)...
    (x>1e-12).*(expint(max(x,1e-12)/tau_2)-expint(max(x,1e-12)/tau_1))+...
    (x<=1e-12).*log(tau_2/tau_1);

intgnd=@(t,tau,tau_1,tau_2,B,C,gamma)...
    (1+C*(delta_E1(t-tau,tau_1,tau_2))).*exp(B*gamma*tau);

%%define eqn5 and eqn6 (as defined in Reference [1])%%
% eqn5 describes the ramping region, eqn6 describes the stress relaxation
% A and B describe the elastic response while C, tau1, and tau2 describe
% the relaxation or viscoelastic response. 
% A: elastic stress constant 
% B: elastic power constant
% C: magnitude of the viscous response (relaxation index)s
% tau1: short relaxation time constant
% tau2: long relaxation time constant
eqn5=@(A,B,C,tau_1,tau_2) arrayfun(@(t_step)...
    ((A*B*gamma)/(1+C*log(tau_2/tau_1)))*...
    integral(@(tau) intgnd(t_step,tau,tau_1,tau_2,B,C,gamma),0,t_step),t_ra);
R=N_ramp/CSA;

eqn6=@(A,B,C,tau_1,tau_2) arrayfun(@(t_step)...
    ((A*B*gamma)/(1+C*log(tau_2/tau_1)))*...
    integral(@(tau) intgnd(t_step,tau,tau_1,tau_2,B,C,gamma),0,t_0),t_rx);
S=N_relax/CSA;

%%nonlinear least-squares regression%%
%levenberg-marquardt was the model used by Abramowitch et al.
min_fun=@(x)...
    [R-eqn5(x(1),x(2),x(3),x(4),x(5));...
    S-eqn6(x(1),x(2),x(3),x(4),x(5))];

options=optimoptions('lsqnonlin','Algorithm','levenberg-marquardt',...
    'MaxIterations',4.00e+3,'Display','iter','FunctionTolerance',1e-6);
[vars,resnorm]=lsqnonlin(min_fun,x0,[],[],options);

%plot for verification
figure()
hold on
plot(t_ra,R,'b');
plot(t_rx,S,'r');
plot(t_ra,eqn5(vars(1), vars(2), vars(3), vars(4), vars(5)),'b');
plot(t_rx,eqn6(vars(1), vars(2), vars(3), vars(4), vars(5)),'r');




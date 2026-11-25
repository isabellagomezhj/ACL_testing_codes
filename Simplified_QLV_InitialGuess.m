%% Simplified Quasi-Linear Viscoelastic Theory (QLV) -- Initial Guess Analysis
% This code utilizes the QLV method found in reference [1] to estimate the
% material properties of a ligament. The QLV method is run for 100 initial
% guesses to determine an optimal initial guess for use in
% "Simplified_QLV_Model.m". Initial guesses are calculated by multiplying
% random numbers spanning 0.1 - 10 by a set of initial coefficients, a
% method consistent with reference [1].
%
% References : 
% [1] "S:\BiomechanicsResearch\groupImhauser\ACL Structure Mechanics\Literature\Abramowitch FiniteStrainRate StressRelaxation 2004.pdf"
%
% This code works specifically for files formatted the same as:
% "S:\BiomechanicsResearch\groupImhauser\ACL Structure Mechanics\Data_Reduced\pcl_07012026.xlsx"
% and will not work for another file that is formatted differently

close all
clear 
clc

fprintf('Trial started %s ',string(datetime('now'))); %print start date/time
run_started=datetime('now'); %record start date/time

% set location of raw data
data="S:\BiomechanicsResearch\groupImhauser\ACL Structure Mechanics\Data_Reduced\pcl_11122026.xlsx";
sheet=sheetnames(data);

% define initial guesses and relaxation duration
x_start=[0.5 50 0.075 0.5 75]; %define starting set of initial guesses
rand_no=0.1+(10-0.1)*rand(1,100)'; %makes list of 100 rand #s 0.1-10

x_guesses=zeros(100,5); %initialize initial guess array
for i=(1:length(rand_no))
    x_guesses(i,:)=rand_no(i)*x_start; %multiply x_start by each rand #
end

% must also uncomment line under %%ACU%% below if you want relax time cut
% (this line is just used to organize the outputs into folders)
rx_duration=14; %duration (mins) set to approx relax length (5,10,14) 

testno=1:length(rand_no); %make array to label each initial guess

% define save locations for results (csv file) and figures (.fig files)
run_folder=sprintf('Test_x0[%.2g %.2g %.2g %.2g %.2g]',x_start);
time_folder=sprintf('%g-min',rx_duration);
save_Dir=fullfile("S:\BiomechanicsResearch\groupImhauser\ACL Structure Mechanics\Data_Reduced\QLV Initial Guess Testing",run_folder,time_folder); %folder containing Figures and Results.csv
save_file=fullfile(save_Dir,"Results.csv"); % file name & .csv ext
fig_Dir=fullfile(save_Dir,"Figures"); %folder name w/in Dir

% ensure past data is not overwritten
if exist(save_Dir,'dir')
    prompt=sprintf('The folder %s already exists. Results will be appended and figures will be overwritten. Proceed?',run_folder);
    warn=questdlg(prompt,'Warning','Yes','No','No');

    if strcmpi(warn,'No')||isempty(warn)
        disp('Please move or rename the existing folder, or change x0');
        return;
    end
end

if ~exist(fig_Dir,'dir'),mkdir(fig_Dir); end %make folders if ~exist

%% Geometric Data
% cross sectional area should be calculated using "CalculatingCSA.m" with
% .stl files collected using the 3D Scanner
% a factor of 1/3 was used to reflect the area of the posteromedial bundle
CSA=56.8236953460548*(1/3);

%% Acquire relevant data from files
sheet_1=readtable(data,'Sheet',sheet(1)); %DIC data (raw)
sheet_2=readtable(data,'Sheet',sheet(2)); %Acumen data (raw)

%% Establish relevant time ranges
% The time ranges for the DIC and Acumen have to be filtered to reflect
% just the ramping and relaxation portion of the testing. Both are sampled
% using a combination of masking and detecting changes from either the raw
% or smoothed data.

%%Acumen%%
% the Acumen data must be cut to include just before the ramp begins
% (before the dip that starts the ramp) and to just after the stress
% relaxation ends (include the decrease that occurs at the end of the
% relaxation plateau). Ensure the time vector for the cut data is adjusted
% to begin at t=0.

ACU_t_raw=sheet_2.TimeZeroed; %sheet w/ filtered and zeroed Acumen time
N_raw=sheet_2.N; %sheet w/ filtered Acumen force outputs

% % plot to view raw data
% figure()
% plot(ACU_t_raw,N_raw)
% title("Raw Data - ACU");

% filter raw data
start=40; %a time after cycles but before ramp to 4%
stop=980; %a time after stress relaxation ends, before ramp to 5%

mask1=(ACU_t_raw>=start)&(ACU_t_raw<=stop);
ACU_t=ACU_t_raw(mask1);
ACU_t=ACU_t-ACU_t(1);
N=N_raw(mask1);

%find the start of the ramp and the end of the relaxation
N_smooth=smoothdata(N,'movmean',25);
TF=find(ischange(N_smooth,'linear'));

[~,k]=mink(N_smooth(TF),2);
ramp_start=min(TF(k)); %min index will be ramp start

relax_end_vals=N_smooth(TF(ACU_t(TF)>920)); %max index will be relax end
relax_end=find(N_smooth==max(relax_end_vals)); 

% find the start of the relaxation
[peaks,locs]=findpeaks(N,ACU_t);
relax_start=find(N==max(peaks));

% define any critical points for the regression
ACU_t_ramp=ACU_t(ramp_start:relax_start);
N_ramp=N(ramp_start:relax_start);

ACU_t_relax=ACU_t(relax_start:relax_end);
N_relax=N(relax_start:relax_end);

t_ra=ACU_t_ramp-ACU_t_ramp(1); %zeroes the ramp time 
t_rx=ACU_t_relax-ACU_t_ramp(1); %matches the relax time to new ramp start
t_0=t_rx(1); %defines t_0 as the relax start (matches ramp end)

% % shorten the relaxation (optional) - comment out if not using
% frames=(rx_duration*60)*100; %frames (assuming 100Hz) in duration
% t_rx=t_rx(1:frames); %redefine the ramp

% % plot to verify smoothing and critical point identification
% fig_ACU=figure()
% hold on
% plot(ACU_t,N);
% plot(ACU_t([ramp_start,relax_start,relax_end]),N([ramp_start,relax_start,relax_end]),'y*');
% plot(ACU_t_ramp,N_ramp,'g');
% plot(ACU_t_relax,N_relax,'m');
% title('Filtered Data (W/ Critical Pts.) - ACU');

% downsample relaxation ONLY (5Hz was used in [1])
N_relax_5=downsample(N_relax,20);
ACU_t_relax_5=downsample(ACU_t_relax,20);

% redefine critical pts. that change with downsampling
t_rx=ACU_t_relax_5-ACU_t_ramp(1);
t_0=t_rx(1);

%%DIC%%
% The code automatically cuts the DIC data for the test on 11/12/2026.
% Adjust the bound of mask 2 to reflect the dataset in use. Ensure the 
% time vector for the cut data begins at t=0 by adjusting the time scale.

ystrain_raw=rmmissing(sheet_1.AverageStrainY);
DIC_t_raw=linspace(0,length(ystrain_raw)/25,length(ystrain_raw)); %calculates time based on 25Hz sample rate

% % plot to view raw data
% figure()
% plot(DIC_t_raw,ystrain_raw)
% title("Raw Data - DIC");

% cut the data to one cycle
mask2=(DIC_t_raw<=35); %select a pt. after one complete cycle
DIC_t1=DIC_t_raw(mask2);
DIC_t1=DIC_t1-DIC_t1(1); %zeroes the filtered data
ystrain1=ystrain_raw(mask2);

% % plot to verify adequate filtering
% figure()
% plot(DIC_t1,ystrain1);
% title('Filtered Data - DIC');

% smooth data and detect changes
ystrain1_smooth=smoothdata(ystrain1,'movmean');
TF=ischange(ystrain1_smooth,'linear');
idxs=find(TF);

% % plot to verify smoothing and critical point identification
% figure()
% hold on
% plot(DIC_t1,ystrain1_smooth);
% plot(DIC_t1(TF),ystrain1(TF),'r*');
% title('Data Smoothed (W/ Critical Pts.) - DIC');

% define ramp regions (used to calculate gamma)
DIC_t_ra=DIC_t1(1:idxs(1));
ystrain_ra=ystrain1(1:idxs(1)).*0.01;

% define relax regions
DIC_t_rx=DIC_t1(idxs(1):idxs(2));
ystrain_rx=ystrain1(idxs(1):idxs(2)).*0.01;

% % plot to verify sectioning
% fig_DIC=figure()
% hold on
% plot(DIC_t_ra,ystrain_ra,'g');
% plot(DIC_t_rx,ystrain_rx,'m');
% plot(DIC_t1(TF),ystrain1(TF)*0.01,'y*');
% title('Filtered Data (W/ Critical Pts.) - DIC');

%% Calculate strain rate (gamma)
coeffs=polyfit(DIC_t_ra,ystrain_ra,1); %fit ramp data (from DIC) linearly
gamma_fit=polyval(coeffs,DIC_t_ra);
gamma=coeffs(1); %slope of best fit line

% % plot to view fit
% fig_gamma=figure()
% hold on
% plot(DIC_t_ra,ystrain_ra,'.');
% plot(DIC_t_ra,gamma_fit);
% title('Gamma');

%% Non-Linear Regression (Levenberg-Marquardt)

% define the change in the exponential integrals
%
% delta_E1 avoids integral() rounding to a zero numerator which results in
% Inf or NaN values in the regression by evaluating small values at the
% finite limit of the difference (E1()-E1()) found in the integrand
delta_E1=@(x,tau_1,tau_2)...
    (x>1e-12).*(expint(max(x,1e-12)/tau_2)-expint(max(x,1e-12)/tau_1))+...
    (x<=1e-12).*log(tau_2/tau_1);

% define eqn5 and eqn6 (as defined in Reference [1]
%
% eqn5 describes the ramping region, eqn6 describes the stress relaxation
% A and B describe the elastic response while C, tau1, and tau2 describe
% the relaxation or viscoelastic response. 
% A: elastic stress constant 
% B: elastic power constant
% C: magnitude of the viscous response (relaxation index)s
% tau1: short relaxation time constant
% tau2: long relaxation time constant
tol={'ArrayValued',true,'RelTol',1e-10,'AbsTol',1e-12};

eqn5=@(A,B,C,tau_1,tau_2) ((A*B*gamma)/(1+C*log(tau_2/tau_1)))*...
    integral(@(u) t_ra.*(1+C*delta_E1(t_ra*(1-u),tau_1,tau_2))...
    .*exp(B*gamma*t_ra*u),0,1,tol{:});
R=N_ramp/CSA; %MPa

eqn6=@(A,B,C,tau_1,tau_2) ((A*B*gamma)/(1+C*log(tau_2/tau_1)))*...
    integral(@(tau) (1+C*delta_E1(t_rx-tau,tau_1,tau_2))...
    .*exp(B*gamma*tau),0,t_0,tol{:});
S=N_relax_5/CSA; %MPa

% levenberg-marquardt was the algorithm used by Abramowitch et al. [1]
min_fun=@(x)...
    [R-eqn5(x(1),x(2),x(3),x(4),x(5));...
    S-eqn6(x(1),x(2),x(3),x(4),x(5))];

options=optimoptions('lsqnonlin','Algorithm','levenberg-marquardt',...
    'MaxIterations',4.00e+3,...
    'FiniteDifferenceStepSize',1e-6,...
    'FunctionTolerance',1e-8,...
    'StepTolerance',1e-6,...
    'Display', 'iter');

for i=1:size(x_guesses,1)
    
    x0=x_guesses(i,:); %set the current initial guess
    test=testno(i);

    [vars,resnorm,~,exitflag]=lsqnonlin(min_fun,x0,[],[],options);
    
    results=[vars resnorm];
    
    % plot for verification
    fig_fit=figure('Visible','off');
    hold on
    plot(t_ra,R,'.k');
    plot(t_rx,S,'.k','HandleVisibility','off');
    plot(t_ra,eqn5(vars(1), vars(2), vars(3), vars(4), vars(5)),'.');
    plot(t_rx,eqn6(vars(1), vars(2), vars(3), vars(4), vars(5)),'.');
    title("Curve-Fit");
    legend('Raw Data','Ramp - Fit','Relax - Fit');

    % save results and figures
    row=table(test,x0(1),x0(2),x0(3),x0(4),x0(5),vars(1),vars(2),vars(3),vars(4),vars(5),resnorm,exitflag,...
    'VariableNames',{'trial','A_0','B_0','C_0','tau1_0','tau2_0','A','B','C','tau1','tau2','resnorm','exitflag'});
    writetable(row,save_file,'WriteMode','append','WriteVariableNames',~isfile(save_file));
    savefig(fig_fit,fullfile(fig_Dir,sprintf("trial_%03d",test)));    
    close(fig_fit)

    fprintf('Trial %d finished %s | resnorm %.4g\n',test,string(datetime('now')),resnorm); % print that the trial has completed
    % run_finished=datetime('now'); %record end date/time

end

savefigpng(fig_Dir,0); %save figs (.png format) in figure directory

% % optional figures to save (must uncomment plots above as well)
% savefig(fig_ACU,fullfile(fig_Dir,sprintf("trial_%03d_ACU_critical_pts",test)));
% savefig(fig_DIC,fullfile(fig_Dir,sprintf("trial_03d_DIC_critical_pts",test)));
% savefig(fig_gamma,fullfile(fig_Dir,sprintf("trial_03d_Gamma_fit",test)));
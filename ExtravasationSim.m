% Anthony Davis & Ryan Knight
% ME 5180 Simulation Project
% 26 Apr. 2016

% Define known values that remain the same

kBT = 4.282; % kBT in pN·nm at physiological temperature (37 deg C)
delta_x = 0.1; % Delta x in nm
N = 1e3; % Number of tumor cell bonds to simulate
t_max = 10; % Length of simulation in s
t = 0:0.01:t_max; % Vector of time values in s
mean_RF = 33; % Mean rupture force in pN from R.H. Eibl & M. Benoit
SD_RF = 12; % Std. Dev. of rupture force in pN from R.H. Eibl & M. Benoit
k_off_0 = 0.0355; % Off rate at F = 0 in s^-1
RF_min = 0; % Minimum rupture force value for histogram
RF_max = 100; % Maximum rupture force value for histogram
bin_size = 10; % Bin size (pN)
RF_bins = (RF_min + bin_size/2):bin_size:(RF_max - bin_size/2);

% Initialize vectors and counter variables

F_rup = zeros(1,N); % Initialized rupture force vector
F_dot = zeros(1,N); % Initialized loading rate vector
F = zeros(N,length(t)); % Initialized force vector
p_surv = zeros(N,length(t)); % Initialized survival probability vector
k_off = zeros(N,length(t)); % Initialized off rate vector
P_f = zeros(N,length(t)); % Initialized rupture force probability vector
Rup_Y = 0; % Counter variable for detached tumor cell bonds
Rup_N = 0; % Counter variable for attached tumor cell bonds

% Simulation is run for 10 seconds of force loading on tumor cell bonds. 
% This process is done for 1000 tumor cell bonds.
% This simulation can be done with different values (parametric study).

for i = 1:N % Begin for loop
    
for j = 1:length(t) % Begin nested for loop
    
    % Tumor cell bonds are given randomly determined rupture forces in pN
    
    F_rup(i) = abs(33+12*randn);
    
    % Loading rate for each tumor cell bond is calculated
    
    F_dot(i) = (k_off_0*kBT/delta_x)*exp(F_rup(i)*delta_x/kBT); 
    
    % Force in pN at each time increment is calculated from loading rate
    
    F(i,j) = F_dot(i)*t(j);
    
    % Survival probability distribution is calculated
    
    p_surv(i,j) = exp(((-k_off_0*kBT)/(F_dot(i)*delta_x))*(exp((F(i,j)*delta_x)/kBT)-1)); 
    
    % Off rate at each time increment is calculated
    
    k_off(i,j) = k_off_0*exp((F(i,j)*delta_x)/kBT);
    
    % Rupture force probability distribution is calculated
    
    P_f(i,j) = (k_off(i,j)*p_surv(i,j))/F_dot(i);
    
end % End nested for loop
    
% Check if the bond has ruptured (i.e. if tumor cell bond has detached)

if F(i,length(t)) <= F_rup(i) % Begin nested if statement
    
    Rup_N = Rup_N+1; % Denote that tumor cell bond has remained attached
    
else 
    
   Rup_Y = Rup_Y+1; % Denote that tumor cell bond has detached
    
end % End nested if statement

end % End for loop

% Create histogram showing rupture force distribution

figure (1)
hold on
RF_hist = hist(F_rup,RF_bins);
bar(RF_bins,RF_hist,'FaceColor',[1 1 1]*0.4);
xlabel('Rupture Force (pN)','FontSize',16)
ylabel('# of Tumor Cell Bonds','FontSize',16)
hold off

% Calculate statistics and relevant values

% Highest rupture force in pN is determined

F_rup_H = max(F_rup);

% Lowest rupture force in pN is determined

F_rup_L = min(F_rup);

% Average rupture force in pN is determined

F_rup_A = mean(F_rup);

% Standard deviation rupture force in pN is determined

F_rup_SD = std(F_rup);

% Highest loading rate in pN/s is determined

F_dot_H = max(F_dot);

% Lowest loading rate in pN/s is determined

F_dot_L = min(F_dot);

% Average loading rate in pN/s is determined

F_dot_A = mean(F_dot);

% Standard deviation loading rate in pN/s is determined

F_dot_SD = std(F_dot);

% Average force vector is determined

F_A = mean(F,1);

% Average survival probability vector is determined

p_surv_A = mean(p_surv,1);

% Average off rate vector is determined

k_off_A = mean(k_off,1);

% Average rupture force probability is determined

P_f_A = mean(P_f,1);

% Show Results of Simulation

% Display survival rate after the time simulated

fprintf('\n\nAfter %d s, %0.1f %% of the tumor cell bonds will still be attached to the endothelium.\n\n',t_max,(Rup_N/N)*100)

% Display average rupture force and average loading rate

fprintf('The average rupture force is %0.1f pN and the average loading rate is %0.2f pN/s.\n\n',F_rup_A,F_dot_A)

% Create plot of survival probabilty as a function of time

figure(2), plot(t,p_surv_A), title('Survival Probability vs. Time')
xlabel('Time (s)'), ylabel('Survival Probability')

% Create plot of survival probabilty as a function of force

figure(3), plot(F_A,p_surv_A), title('Survival Probability vs. Force')
xlabel('Force (pN)'), ylabel('Survival Probability')

% Create plot of left half of rupture force probability distribution

figure(4), plot(F_A,P_f_A), title('Rupture Force Probabilty Distribution - Left Half of Curve')
xlabel('Force (pN)'), ylabel('Rupture Force Probability')


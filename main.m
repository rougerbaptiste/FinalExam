close all
clear all 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Analysis of the original model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_names= {'N','E','A'}; % N is in CFU
p_names= {'k','Nm','d','ke','de','va','da'};

%define reference parameter values (from paper)
p_ref([3 4 5 6])= [4e-03 5 2 4.8e-07]; %set values for d, ke, de, and va
p_ref([1 2 7])= [0.970 1.24e09 0.639];% set values for k, Nm and da at pH 7

%define initial condifions and simulation time (somewhat arbitrary)
x0= [100 0 0];
tspan= [0:1:60];

%In order to avoid doing always the same computations, you can set several Q to 0...
Q1=1 ; Q2=1; Q3=1; Q4=1; Q5=1; Q6_7=1; Q8=1; Q9=0; Q10=0;

% Computes ON and OFF behaviors of engineered system; plot on same figure
% %
if Q1
    display('--------------Question 1--------------');
    %set parameters p for the ON system and do numerical simulation
    p= p_ref;
    [t x_on]= ode15s(@you_ode,tspan, x0, [], p);
    %set parameters p for the OFF system and do numerical simulation
    p([4 6]) = [0 0]; %define p
    [t x_off]= ode15s(@you_ode,tspan, x0, [], p);
    figure(1)
    subplot(2,1,1); plot(t,x_on(:,1)/p_ref(2)); legend('N/Nm'); hold on 
    subplot(2,1,2); plot(t,x_on(:,2:3)); legend('E','A'); hold on
    subplot(2,1,1); plot(t,x_off(:,1)/p_ref(2),'--'); legend('N/Nm');
    subplot(2,1,2); plot(t,x_off(:,2:3),'--'); legend('E','A');
end

% % %
% % Compute cell density of engineered system at pH6.2 and pH7.8 and the corresponding fold change
% % %

p_pH62= [0.885 1.25e09 4e-03 5 2 4e-07 0.274];% stores values of k, Nm and da at pH 6.2   
p_pH78= [0.936 1.20e09 4e-03 5 2 4e-07 1.19];% same at pH 7.8   
if Q2
    display('--------------Question 2--------------');
    %compute cell density at steady state for pH6.2; stored in Nmin
    p = p_pH62; %define p
    [t x]= ode15s(@you_ode,tspan, x0, [], p);
    Nmin= x(end,1) %define Nmin using x
    %compute cell density at steady state for pH7.8; stored in Nmax
    p = p_pH78; %define p
    [t x]= ode15s(@you_ode,tspan, x0, [], p);
    Nmax= x(end,1) %define Nmax using x
    %define and display the fold change associated to pH variations
    fold_change_pH= Nmax/Nmin % define fold_change_pH as a function of Nmin and Nmax;
    display(['fold change; ' num2str(fold_change_pH)]); % just for nice display; "num to string" converts numbers to strings of characters
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Analysis of extended  models
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %extended parameter set
% %all parameters unchanged, excepted m, theta, and eta  
% pe_names= {'k','Nm','d','ke','de','va','da','m','theta','eta'};

if Q3
    display('--------------Question 3: see text--------------');
end
if Q4
    display('--------------Question 4: see you_ode, you_odeR, you_odeI, and you_odeRI functions--------------');
end

%pe stands for parameters of extended system
pe_ref(1:7)= p_ref;
pe_ref([8 9 10])= [1 1 2.1]; %m=1, theta=1, and eta= 2.1

if Q5
    display('--------------Question 5--------------');
    figure(2)
    pe= pe_ref; %define pe
    [t x]= ode15s(@you_odeR,tspan, x0, [], pe); % compute behavior of circuit R
    subplot(2,1,1); plot(t,x(:,1)/pe_ref(2)); legend('N/Nm');
    subplot(2,1,2); plot(t,x(:,2:3)); legend('E','A');
    ss= x(end,:); % store steady state values for N, E, and A
    display(['ss for circuit R: ' num2str(ss)]);
    figure(3)
    pe= pe_ref;
    [t x]= ode15s(@you_odeI,tspan, x0, [], pe); % compute behavior of circuit I
    subplot(2,1,1); plot(t,x(:,1)/pe_ref(2)); legend('N/Nm');
    subplot(2,1,2); plot(t,x(:,2:3)); legend('E','A');
    ss= x(end,:);
    display(['ss for circuit I: ' num2str(ss)]);
    figure(4)
    pe= pe_ref;
    [t x]= ode15s(@you_odeRI,tspan, x0, [], pe); % compute behavior of circuit RI
    subplot(2,1,1); plot(t,x(:,1)/pe_ref(2)); legend('N/Nm');
    subplot(2,1,2); plot(t,x(:,2:3)); legend('E','A');
    ss= x(end,:);
    display(['ss for circuit RI: ' num2str(ss)]);
end

m= [0 0.1 0.2 0.4 0.6 0.8 1 1.5 2 3];
%m= [0 1 3]; %for debugging purpose, you can do computations with this m vector. computations are fastest!
if Q6_7
    display('--------------Question 6 and 7--------------');
    %compute nominal I/O behavior for the 3 models
    pe= pe_ref; % notably theta=1 and eta= 2.1. the value of m will be changed in computeIO_X
    figure(5)
    Nnom= computeIO_N(m, tspan, x0, pe);
    semilogy(m, Nnom); legend('N you', 'N you R', 'N you I', 'N you RI'); hold on
    figure(6)
    Enom= computeIO_E(m, tspan, x0, pe);
    semilogy(m, Enom); legend('E you', 'E you R', 'E you I', 'E you RI'); hold on
    figure(7)
    Anom= computeIO_A(m, tspan, x0, pe);
    semilogy(m, Anom); legend('A you', 'A you R', 'A you I', 'A you RI'); hold on

    %compute mean square deviation with respect to nominal behavior for the 3 models and for specific parameter values
    pe= pe_ref;
    pe([9 10])= [0.4 2.4];
    %msd(1) stores the msd for E in model R
    %msd(2) stores the msd for A in model R
    %msd(3) stores the msd for E in model I
    %msd(4) stores the msd for A in model I
    %msd(5) stores the msd for E in model RI
    %msd(6) stores the msd for A in model RI
    theta = 0.4;
    eta = 2.4;
    msd(1)= compute_msd(m,tspan,x0,pe,pe_ref,1,1);
    msd(2)= compute_msd(m,tspan,x0,pe,pe_ref,1,2);
    msd(3)= compute_msd(m,tspan,x0,pe,pe_ref,2,1);
    msd(4)= compute_msd(m,tspan,x0,pe,pe_ref,2,2);
    msd(5)= compute_msd(m,tspan,x0,pe,pe_ref,3,1);
    msd(6)= compute_msd(m,tspan,x0,pe,pe_ref,3,2);
    display(['msd: ' num2str(msd)]);
    %display(['msd : ' msd1 msd2 msd3 msd4 msd5 msd6]);
end

if Q8
    display('--------------Question 8--------------');
    %theta= [0.4 0.6 0.8 1 1.2 1.4 1.6];
    %eta= [1.8 2 2.2 2.4];
    theta= [0.4 1 1.6]; % again for debugging, you can use these values instead
    eta= [1.6 2 2.4];
    msd= zeros(length(theta),length(eta));
    
    for i= 1:length(theta)
        for j= 1:length(eta)
            pe= pe_ref;
            pe([9 10])= [theta(i) eta(j)];
            msd(i,j)= compute_msd(m,tspan,x0,pe,pe_ref,1,1); % here store in msd(i,j) the mean square deviation between the nominal IO behavior and the one one obtains with specific values of theta and eta, if the model is youR, and the observed variable is E
            %use compute_msd and pe
        end
    end
    variance_IO(1)= var(msd(:));

    for i= 1:length(theta)
        for j= 1:length(eta)
            pe= pe_ref;
            pe([9 10])= [theta(i) eta(j)];
            msd(i,j)= compute_msd(m,tspan,x0,pe,pe_ref,1,2); % the same as above, but if the model is youR and the observed variable is A
        end
    end
    variance_IO(2)= var(msd(:));
    
    for i= 1:length(theta)
        for j= 1:length(eta)
            pe= pe_ref;
            pe([9 10])= [theta(i) eta(j)];
            msd(i,j)= compute_msd(m,tspan,x0,pe,pe_ref,2,1); % the same as above, but if the model is youI and the observed variable is E
        end
    end
    variance_IO(3)= var(msd(:));
    
    for i= 1:length(theta)
        for j= 1:length(eta)
            pe= pe_ref;
            pe([9 10])= [theta(i) eta(j)];
            msd(i,j)= compute_msd(m,tspan,x0,pe,pe_ref,2,2); % the same as above, but if the model is youI and the observed variable is A
        end
    end
    variance_IO(4)= var(msd(:));
    
    for i= 1:length(theta)
        for j= 1:length(eta)
            pe= pe_ref;
            pe([9 10])= [theta(i) eta(j)];
            msd(i,j)= compute_msd(m,tspan,x0,pe,pe_ref,3,1); % the same as above, but if the model is youRI and the observed variable is E
        end
    end
    variance_IO(5)= var(msd(:));
    
    for i= 1:length(theta)
        for j= 1:length(eta)
            pe= pe_ref;
            pe([9 10])= [theta(i) eta(j)];
            msd(i,j)= compute_msd(m,tspan,x0,pe,pe_ref,3,2); % the same as above, but if the model is youRI and the observed variable is A
        end
    end
    variance_IO(6)= var(msd(:));
    
    display(['variance: ' num2str(variance_IO)]); 
    % make you choice here!
end

% % experimental data should be stored in data.mat file. This file must be in the same folder as your main file
% if Q9



%     load('data.mat','data');
%     figure(8)
%     plot(..,..,'x'); legend('E or A'); hold on  %here, plot the experimental data 
%     k(:,1)= pe_ref([.. ..]); %here data for reference values for theta and eta
%     k_sigma= k/3;
%     opts=cmaes;
%     opts.DispModulo= 20; %displays info every 20 iterations (default is 100)
%     opts.StopFitness= 1e-2; %displays info every 20 iterations (default is 100)
%     [k_opt, cost_min, counteval, stopflag, out, bestever]= cmaes('compute_cost', k, k_sigma,opts,m, tspan, x0, pe_ref, data);
%     
%     ..%here define pe such that pe stores the parameters you found
%     ..
%     for i=1:length(m)
%         .. %here pe stores the ith inducer concentration
%         [t x]= ode15s(..);
%         .. %store in A(i) or E(i) the steady state values for A or E that you computed
%     end
%     .. % plot A or E (the plot goes in figure 8)
%     save('k_opt.mat','k_opt'); %save these values to avoid computing them each time
% end
% if Q10
%     load('k_opt.mat', 'k_opt'); %reload the optimal parameters found
%     .. % define pe_opt based on pe_ref and k_opt
%     ..
%     .. % set m=0 in pe_opt  -- line A (for future use)
%     [t x]= ode15s(..); %do the num sim for the you model
%     max(1)= .. %stores the cell density at steady state in max(1)
%     .. % set m=3 in pe_opt
%     [t x]= ode15s(..); %do the num sim for the you model
%     min(1)= ..%stores the cell density at steady state in min(1) -- line B (for future use)
%     .. %redo from line A to line B with the youR model, storing info in max(2) and min(2)
%     .. %redo from line A to line B with the youI model, storing info in max(3) and min(3)
%     .. %redo from line A to line B with the youRI model, storing info in max(4) and min(4)
%     figure(9)
%     plot(1:4, max./min,'o');
% 
%     .. % compute timed behavior of best model in absence of inducer (behavior stored in x1) or in presence of inducer (behavior stored in x2) 
%     ..
%     ..
%     ..
%     figure(10)
%     semilogy(t,x1(:,1)/p_ref(2),'b'); legend('N/Nm'); hold on
%     semilogy(t,x2(:,1)/p_ref(2),'k'); legend('N/Nm'); hold on
% end
   
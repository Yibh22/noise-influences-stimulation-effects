function [signals] = wc_coupled_stochastic1_sd5(G,D,time,dt,c5,stim_P,sigma,res,rs)



%% Equation and simulation parameters: 
N = size(G,2);                      % --- Number of brain areas
noise_var = sigma;                % --- Noise amplitude



initEstate = 0.0*randn(1,N)+0.1;    % --- size 1xN, gives the inital states of the 
                                          %excitatory population for each oscillator
initIstate = 0.0*randn(1,N)+0.1;    % --- size 1xN, gives the initial states of the
                                          %inhibitory population for each oscillator

% Model parameters 
tau = 8;            % --- the excitatory time constant in ms

c1 =  16;           % --- parameter that defines the self-coupling of the
                    %     excitatory population
c2 = 12;            % --- parameter that defines the coupling of the
                          %inhibitory to excitatory population for a single 
                          %column
c3 = 15;            % --- parameter that defines the coupling of the
                          %excitatory to inhibitory population for a single 
                          %column
c4 = 3;             % --- parameter that defines the self-coupling of the
                          %inhibitory population
a_e = 1.3;          % --- parameter for excitatory sigmoidal
a_i = 2;            % --- parameter for inhibitory sigmoidal
theta_e = 4;        % --- parameter for excitatory sigmoidal
theta_i = 3.7;      % --- parameter for inhibitory sigmoidal


%% Setting up delay matrix, converting time delay matrix into index
delays_m_steps=round(D./dt); 


%% Simulation parameters

%setting the self_connectivity to 0 
for i=1:N
    G(i,i)=0;
end

% total simulation steps
T = (time/dt)+1;
signals.t = transpose(0:dt:time);

%parameters for sigmoid function
Smax_e=1-1/(1+exp(a_e*theta_e));
Smax_i=1-1/(1+exp(a_i*theta_i));
sig_params_e=[a_e theta_e];
sig_params_i=[a_i theta_i];


% initializing output (setting all the initial values same) 

estate = zeros(T, N);
istate= zeros(T, N);

estate(1,:) = initEstate;
istate(1,:) = initIstate;



%% Setting up stimulation matrix

%Excitatory stimulatioin (P)
P = zeros(1,N);          
P_matrix=repmat(P,T,1);  %T*N
num_stim=size(stim_P);
num_stim=num_stim(1); %n
if num_stim > 0
    for i=1:num_stim
        oscillator_i=stim_P(i,1);
        amp_stim_i=stim_P(i,2);
        start_stim_i=stim_P(i,3)/dt;
        end_stim_i=stim_P(i,4)/dt;
        P_matrix(start_stim_i:end_stim_i,oscillator_i)=P_matrix(start_stim_i:end_stim_i,oscillator_i)+amp_stim_i;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main Loop: Simulating coupled and stimulated WC oscillators 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for t= 2:T
    %creating noise that is fed to the excitatory and inhibitory populations
    noise_e=noise_var * randn(rs,1,N);
    noise_i=noise_var * randn(rs,1,N);
    
    %calculating the delay time
    delayed_time_index=t-delays_m_steps;
    negative_index=find(delayed_time_index < 2);
    if ~isempty(negative_index)
        delayed_time_index(negative_index)=2;
    end
    
    %calcuating the estate and istate of the delay times
    edelayed=zeros(N);
    for i=1:N
        edelayed_i=estate(delayed_time_index(:,i)-1,i);
        edelayed(i,:)=edelayed_i;
        
    end
   
    
    %calculating the coupling term input
    c_in_e = sum(G.*edelayed,1);

    
    
       
    %% Calculating the states using Euler-Maruyama method 

    estate(t,:) = estate(t-1,:) + dt .* (1/tau) .* (- estate(t-1, :) + ...
                (Smax_e-estate(t-1,:)).*sigmoidal_wc72(c1.*estate(t-1,:)-c2.*istate(t-1,:)+c5.*c_in_e+P_matrix(t,:), sig_params_e))...
                + ((1/tau).*dt^0.5.*noise_e);

    istate(t,:) = istate(t-1,:) + dt .* (1/tau) .* (- istate(t-1, :)+ ...
                (Smax_i-istate(t-1,:)).*sigmoidal_wc72(c3.*estate(t-1,:)-c4.*istate(t-1,:), sig_params_i))... %%+c6.*c_in_i
                + ((1/tau).*dt^0.5.*noise_i);

end
E_cutd = WC_downsampling(estate',res); 
I_cutd = WC_downsampling(istate',res); 
signals.e = E_cutd; 
signals.i = I_cutd;




function [sig_val] = sigmoidal_wc72( x, sig_params )
%this is to calculate the sigmoidal funciton used in the Wilson-Cowan
%model described in Wilson and Cowan's 1972 paper

%Input:
%   'x'             --- 1xN vector where N is the
%                       number of oscillators
%   'sig_params'    --- 1x2 vector listing the 'a','theta' parameter
%                       values as in the paper

%Output: 
%   'sig_val'       --- the value of the function

a=sig_params(1);
theta=sig_params(2);

sig_val = 1./(1 + exp(- a.*(x - theta))) - 1/(1 + exp(a*theta));
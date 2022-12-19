clear all;
clc;
N_core = 41; % user define,core of computer
save_dir_input=pwd; %= <user define their own output directory>;
main_dir = pwd; 
data_dir1 = fullfile(main_dir,'Aij_human_Streamlines_scale82.mat'); %connect the path
data_dir2 = fullfile(main_dir,'Dij_human_Streamlines_scale82.mat'); %connect the path
save_dir = save_dir_input;
cluster_dir = save_dir_input;
cluster = parcluster('local'); 
cluster.JobStorageLocation = cluster_dir;
parpool(cluster,N_core);
load(data_dir1)
load(data_dir2)
G = Aij;
D = Dij/10*2.4;
time = 3000; %ms
dt = 0.005;
res = 200; %resolution
c5 = 0.1;
c6 = 0.1;
fix = 1:9;
noise1 = fix*0.000000001;
noise2 = fix*0.00000001;
noise3 = fix*0.0000001;
noise4 = fix*0.000001;
noise5 = fix*0.00001;
noise6 = fix*0.0001;
noise7 = fix*0.001;
noise8 = 0.01;
noise = [noise1,noise2,noise3,noise4,noise5,noise6,noise7,noise8];
s = RandStream.create('mrg32k3a','NumStreams',157440,'Seed','shuffle','CellOutput',true);
for k = 1:30
    for j = 1:64
        sigma = noise(j);
        result = {};
        parfor i = 1:82  
            tic
            rnum = (k-1)*64*82+(j-1)*82+i;
            rs = s{rnum};
            stim_P = [i,1.25,2001,3000];
            [signals] = wc_coupled_stochastic1_sd5(G,D,time,dt,c5,stim_P,sigma,res,rs);
            result{i,1} = signals.e;
            disp(j)
            disp(i)
            toc
        end
        parsave_Euler(save_dir,result,j,k)
    end
end
poolobj = gcp('nocreate');
delete(poolobj);
function [] = parsave_Euler(dir,x,j,k)
%save x,y in dir
% so I can save in parfor loop
d = [dir '/1nodePen20bn_node',num2str(k),'noise',num2str(j),'_norepeat','s0005_Euler.mat'];

save(d,'x');

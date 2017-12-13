clear

addpath(genpath(pwd));
input = 'list.txt';

theta = 0.2;	

ih = fopen(input,'r');

while ~feof(ih)
    file_abspath = fgets(ih); % read line by line
    file_abspath = file_abspath(1:(end-1)); % Final \n probably caused a problem
    fprintf('%s\n',file_abspath);
    
    [~,file_name,~] = fileparts(file_abspath);
    outdir = sprintf('%s_results',file_name);
    system(sprintf('mkdir -p %s',outdir));
    N = preprocessing(file_abspath,outdir,theta);
    decimation(N,sprintf('%s/msa_numerical.txt',outdir),sprintf('%s/weights.txt',outdir),sprintf('%s/decimation_results',outdir));
end
fclose(ih);

clear

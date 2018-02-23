function decimation(N,input_alignment, input_weights, outdir)

	q = 21;
	Nit = 41; % Max. number of decimation rounds
	r = 0.025; % Fraction (approximate) of decimated couplings at each round.
	lambda = 0.01; % Regularization for plm.
	theta = 0.2; % Reweighting threshold. -- Unused here
	ncores = 1; % Number of cores each plm run should use. 
	
	% Creating out directory
	system(sprintf('mkdir -p %s',outdir));

	% Mapping from [i j] to list of pairs
	pos = 1;
	for i = 1:N
		for j = (i+1):N
			pairlist(pos,:) = [i j];
			pos = pos +1;
		end
	end

	%% Initial inference
	decim = 0;
	mask_file = sprintf('%s/mask0.txt',outdir);
	mask = [];
	dlmwrite(mask_file,mask,'delimiter',' ');
	scores_file = sprintf('%s/score0.txt',outdir);
	parameters_file = sprintf('%s/plmInf_0_mat.txt',outdir);
	plmDCA_asymmetric_mask(input_alignment,scores_file,parameters_file,mask_file,theta,ncores,lambda,input_weights);
	
	
	
	%% Decimation
	for s = 1:Nit
		fprintf('s = %d out of %d\n',s,Nit);
		% Choosing decimated pairs
		Fc = dlmread(scores_file);
		[sFc, rFc] = sort(Fc(:,3),'ascend');
		decim = min(size(pairlist,1), decim + round(r*numel(rFc)));
		fprintf('decim = %d -- numpairs = %d\n', decim, size(pairlist,1));
		mask = [];
		for d = 1:decim 
			mask = [mask ; pairlist(rFc(d),:)];
		end
		mask_file = sprintf('%s/mask%d.txt',outdir,s);
		dlmwrite(mask_file,mask,'delimiter',' ');
	
		% Inference
		scores_file = sprintf('%s/score%d.txt',outdir,s);
		parameters_file = sprintf('%s/plmInf_%d_mat.txt',outdir,s);
		plmDCA_asymmetric_mask(input_alignment,scores_file,parameters_file,mask_file,theta,ncores,lambda,input_weights);
	end

	%% Site-specific decimation
	system(sprintf('./sources/script_mask.sh %s/score0.txt %d %d',outdir, N, floor(N/2)));
	system(sprintf('mv mask_t20.txt %s ; mv mask_tLh.txt %s',outdir,outdir));
	% t20
	scores_file = sprintf('%s/score_t20.txt',outdir);
	parameters_file = sprintf('%s/plmInf_t20_mat.txt',outdir);
	mask_file = 'mask_t20.txt';
	plmDCA_asymmetric_mask(input_alignment,scores_file,parameters_file,mask_file,theta,ncores,lambda,input_weights);
	% tL2
	scores_file = sprintf('%s/score_tLh.txt',outdir);
	parameters_file = sprintf('%s/plmInf_tLh_mat.txt',outdir);
	mask_file = 'mask_tLh.txt';
	plmDCA_asymmetric_mask(input_alignment,scores_file,parameters_file,mask_file,theta,ncores,lambda,input_weights);
end
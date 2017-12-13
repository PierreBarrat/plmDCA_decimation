function decimation(N,input_alignment, input_weights, outdir)

	q = 21;
	Nit = 3; % Max. number of decimation rounds
	r = 0.3; % Fraction (approximate) of decimated couplings at each round.
	lambda = 0.01; % Regularization for plm.
	theta = 0.2; % Reweighting threshold. -- Unused here
	ncores = 3; % Number of cores each plm run should use. 
	
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
end
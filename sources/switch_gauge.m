function  [Jo,ho] = switch_gauge(Ji,hi,mode,q,wt01)


	if strcmpi(mode,'LatticeGas')
		N = max(size(hi))/q;
		ho = zeros(1,q*N);
		Jo = zeros(q*N,q*N);

		Jo = Ji;
		ho = hi;
		fprintf('Hello\n')
		N
		for i = 1:N
			for j = 1:N
				Jo((i-1)*q+(1:q),(j-1)*q+(1:q)) = Jo((i-1)*q+(1:q),(j-1)*q+(1:q)) - repmat(Ji((i-1)*q+(1:q),(j-1)*q+q),1,q);
				Jo((i-1)*q+(1:q),(j-1)*q+(1:q)) = Jo((i-1)*q+(1:q),(j-1)*q+(1:q)) - repmat(Ji((i-1)*q+q,(j-1)*q+(1:q)),q,1);
				Jo((i-1)*q+(1:q),(j-1)*q+(1:q)) = Jo((i-1)*q+(1:q),(j-1)*q+(1:q)) + repmat(Ji((i-1)*q+q,(j-1)*q+q),q,q);	
				ho((i-1)*q+(1:q)) = ho((i-1)*q+(1:q)) + Ji((i-1)*q+(1:q),(j-1)*q+q)';
				ho((i-1)*q+(1:q)) = ho((i-1)*q+(1:q)) - repmat(Ji((i-1)*q+q,(j-1)*q+q),1,q);
			end
        end

		for i = 1:N
			ho((i-1)*q+(1:q)) = ho((i-1)*q+(1:q)) - repmat(hi((i-1)*q+q),1,q);% repmat(ho((i-1)*q+q),1,q);
		end

		% ho(1)=1000;
	elseif strcmpi(mode,'regl2')
		N = max(size(hi))/q;
		ho = zeros(1,q*N);
		Jo = zeros(q*N,q*N);

		[Ji,hi] = switch_gauge(Ji,hi,'0sum',q);
		u = hi / (q+N-1);

		for i = 1:N
			ho((i-1)*q+(1:q)) = hi((i-1)*q+(1:q)) - (N-1)*u((i-1)*q+(1:q));
			for j = (i+1):N
				Jo((i-1)*q+(1:q),(j-1)*q+(1:q)) = Ji((i-1)*q+(1:q),(j-1)*q+(1:q)) + repmat(u((i-1)*q+(1:q))',1,q) + repmat(u((j-1)*q+(1:q)),q,1);
			end
		end
		Jo = Jo + Jo';


	elseif strcmpi(mode,'0sum')
	
		% Works for q = 2 at least -- and q = 21

		N = max(size(hi))/q;
		ho = zeros(1,q*N);
		Jo = zeros(q*N,q*N);

		for i = 1:N
			for j = (i+1):N
				Jo((i-1)*q+(1:q),(j-1)*q+(1:q)) = Ji((i-1)*q+(1:q),(j-1)*q+(1:q)) - repmat(mean(Ji((i-1)*q+(1:q),(j-1)*q+(1:q)),1),q,1) - repmat(mean(Ji((i-1)*q+(1:q),(j-1)*q+(1:q)),2),1,q) + mean(mean(Ji((i-1)*q+(1:q),(j-1)*q+(1:q))));
				ho((i-1)*q+(1:q)) = ho((i-1)*q+(1:q)) + mean(Ji((i-1)*q+(1:q),(j-1)*q+(1:q)),2)' - mean(mean(Ji((i-1)*q+(1:q),(j-1)*q+(1:q))));
				% ho((i-1)*q+(1:q)) = hi((i-1)*q+(1:q)) + mean(Ji((i-1)*q+(1:q),(j-1)*q+(1:q)),2)' - mean(mean(Ji((i-1)*q+(1:q),(j-1)*q+(1:q))));
				% ho((j-1)*q+(1:q)) = hi((i-1)*q+(1:q)) + mean(Ji((i-1)*q+(1:q),(j-1)*q+(1:q)),1);
			end
		end
		Jo = Jo + Jo';

		for i = 1:N
			ho((i-1)*q+(1:q)) = hi((i-1)*q+(1:q));
			% for j = 1:N
			% 	ho((i-1)*q+(1:q)) = ho((i-1)*q+(1:q)) + mean(Ji((i-1)*q+(1:q),(j-1)*q+(1:q)),2)' - mean(mean(Ji((i-1)*q+(1:q),(j-1)*q+(1:q))));
			% end
			ho((i-1)*q+(1:q)) = ho((i-1)*q+(1:q)) - mean(hi((i-1)*q+(1:q)));
		end

	elseif strcmpi(mode,'wt')
		N = max(size(hi))/q;

		Jo = Ji;
		ho = hi;
		for i = 1:N
			alpha = find(wt01((i-1)*q+(1:q))==1);
			for j = 1:N
				beta = find(wt01((j-1)*q+(1:q))==1);
				Jo((i-1)*q+(1:q),(j-1)*q+(1:q)) = Jo((i-1)*q+(1:q),(j-1)*q+(1:q)) - repmat(Ji((i-1)*q+(1:q),(j-1)*q+beta),1,q);
				Jo((i-1)*q+(1:q),(j-1)*q+(1:q)) = Jo((i-1)*q+(1:q),(j-1)*q+(1:q)) - repmat(Ji((i-1)*q+alpha,(j-1)*q+(1:q)),q,1);
				Jo((i-1)*q+(1:q),(j-1)*q+(1:q)) = Jo((i-1)*q+(1:q),(j-1)*q+(1:q)) + repmat(Ji((i-1)*q+alpha,(j-1)*q+beta),q,q);	
				ho((i-1)*q+(1:q)) = ho((i-1)*q+(1:q)) + Ji((i-1)*q+(1:q),(j-1)*q+beta)';
				ho((i-1)*q+(1:q)) = ho((i-1)*q+(1:q)) - repmat(Ji((i-1)*q+alpha,(j-1)*q+beta),1,q);
			end
        end

		for i = 1:N
			alpha = find(wt01((i-1)*q+(1:q))==1);
			ho((i-1)*q+(1:q)) = ho((i-1)*q+(1:q)) - repmat(hi((i-1)*q+alpha),1,q);% repmat(ho((i-1)*q+q),1,q);
		end
	end


end
%   Copyright (C) 2022  Ma Shisong Lab  All Rights Reserved.

classdef rgsGGM
	properties
		gene_num = [];
		gene_name = [];
		pcor_all = [];
		pcor_sampling_num = [];
		samples_num = [];
		RoundNumber = [];
		SigEdges = [];
	end

	methods (Static)
		function obj=rgsGGM (x, round_num, gene_name)

			selected_num = 2000;
			cut_off = 0.02;

			[n,p] = size(x);

			gene_id = 1:p;

			if nargin > 2
				gene_id = regexp( gene_name,'([a-zA-Z0-9_\.-]+)','once','match');
			end
			
			obj.samples_num = n;
			obj.gene_num = p;
			obj.RoundNumber = round_num;
			obj.gene_name = gene_id;

			pcor_all = ones( p, p );
			pcor_sampling_num = zeros( p, p );

			time_trend = zeros(100,1);

			fprintf('Calculating pcor in %d iterations.\n', round_num);

			for i = 1 : round_num
				loop_start_t = clock;
				j = datasample(1:p, selected_num, 'Replace', false);
				sx = x(:,j);
				cov_x = cov (sx);
				ix = inv( cov_x );
				d = diag(sqrt(diag(ix)));
				d = inv( d );
				pc = - d * ix * d;
				pc = eye( selected_num ) * 2 + pc;

				t_idx = mod(i,100) + 1;

				for m = 1 : selected_num
					for n = 1: selected_num
						r = j(m);
						s = j(n);
						if r > s
							pcor_sampling_num(r,s) = pcor_sampling_num(r,s) + 1;
							pcor_mn = pc(m,n);
							if abs(pcor_mn) < abs(pcor_all(r,s))
								pcor_all(r,s) = pcor_mn;
							end
						end
					end
				end

				loop_time = etime(clock, loop_start_t);
				idx_time  = mod(i,100) + 1;
				time_trend(idx_time) = loop_time;
				average_loop_time = mean(time_trend);
				time_left = (round_num - i) * average_loop_time / 3600;

				if i == 100
					fprintf('Estimated to complete in %.2f hours.\n', time_left);
				end

				if mod(i,200) == 0 & i < round_num
					fprintf('%d iterations done. %.2f sec/iteration in average. %.2f hours to go.\n', i, average_loop_time, time_left);
				end

				if i == round_num
					fprintf('%d iterations done.\n',i);
				end

			 end

			 obj.pcor_sampling_num = int16(pcor_sampling_num);

			 rho = corr(x);
			 idx = find(pcor_sampling_num == 0);
			 pcor_all(idx) = 0;	
			 idx = find(pcor_all >= cut_off & pcor_all < 1);
			 [e2,e1] = ind2sub(size(pcor_all), idx);
			 e1 = gene_id(e1);
			 e2 = gene_id(e2);
			 e3 = pcor_all(idx);
			 e4 = pcor_sampling_num(idx);
			 e5 = rho(idx);
			 colName = {'GeneA','GeneB','Pcor','SamplingTime','r'};
			 obj.SigEdges = table(e1,e2,e3,e4,e5,'VariableNames',colName);
			 obj.pcor_all = pcor_all;
		 end
	 end
 end



			

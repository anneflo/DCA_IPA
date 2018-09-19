function [DCA_couplings_zeroSum, Meff] = Compute_DCA_model(encoded_focus_alignment)

q=21;
theta = 0.3;
pseudocount_weight = 0.5;

%compute the empirical frequencies from the sequence data
[Meff, Pij_true, Pi_true, alignment_width] = count_alignment(encoded_focus_alignment, theta, q);

%include pseudocounts
[Pij, Pi] = with_pc(Pij_true, Pi_true, pseudocount_weight, alignment_width, q);

%compute the L*(q-1)*L*(q-1) matrix of correlations
C = Compute_C_mat(Pij, Pi, alignment_width, q); 

%invert the matrix of correlations to get the approximate matrix of direct couplings
invC = inv(C);

%now make gauge change and store the direct couplings
DCA_couplings_zeroSum = zeros(alignment_width,alignment_width,q,q);
for i=1:alignment_width
	for j=i:alignment_width % only fill in the upper triangle of W, because only this part is used (residue pairs where i<j)                   
        % get the block of the correlation matrix that corresponds to sites i and j
        DCA_couplings = Compute_DCA_couplings(invC, i, j, q);		
        % change the gauge in this block to the zero-sum gauge
        DCA_couplings_zeroSum(i,j,:,:) = Change_gauge(DCA_couplings,q);        
	end
end

end



%%auxiliary function definitions

function [Meff, Pij_true, Pi_true, alignment_width, q] = count_alignment(encoded_focus_alignment, theta, q)

%calculate Meff
[alignment_height,alignment_width] = size(encoded_focus_alignment);
W = ones(1, alignment_height);
if(theta > 0.0) %whether you weight or not
    W = (1./(1+sum(squareform(pdist(encoded_focus_alignment, 'hamm')<theta))));
end
Meff=sum(W);

%compute the frequencies
Pij_true = zeros(alignment_width, alignment_width, q, q); %a 4-dim matrix (equivalent to a q^2*alignment_width^2 square matrix)
Pi_true = zeros(alignment_width, q);

%single-site frequencies
for j=1:alignment_height %sequence index in the sequence alignment
	for i=1:alignment_width %site index in the sequence alignment
		Pi_true(i, encoded_focus_alignment(j, i)) = Pi_true(i, encoded_focus_alignment(j, i)) + W(j); %increment the proba to have this residue at position i by the weight W(j) of the sequence j considered
	end
end
Pi_true = Pi_true/Meff; %normalization

%two-site frequencies
for l=1:alignment_height %sequence index
	for i=1:alignment_width-1 %site index
		for j=i+1:alignment_width %site index
			Pij_true(i, j, encoded_focus_alignment(l, i), encoded_focus_alignment(l, j)) = Pij_true(i, j, encoded_focus_alignment(l, i), encoded_focus_alignment(l, j)) + W(l);
			Pij_true(j, i, encoded_focus_alignment(l, j), encoded_focus_alignment(l, i)) = Pij_true(i, j, encoded_focus_alignment(l, i), encoded_focus_alignment(l, j));
		end
	end
end
Pij_true = Pij_true/Meff;

%fill the diagonal of Pij_true
scra = eye(q, q); %has to be the same amino acid
for i=1:alignment_width
	for alpha=1:q
		for beta=1:q
			Pij_true(i, i, alpha, beta) = Pi_true(i, alpha) * scra(alpha, beta);
		end
	end
end

end



function [Pij, Pi] = with_pc(Pij_true, Pi_true, pseudocount_weight, alignment_width, q)
%add pseudocounts to deal with some finite size effects 

Pij = (1.-pseudocount_weight)*Pij_true + pseudocount_weight/q/q*ones(alignment_width, alignment_width, q, q);
Pi = (1.-pseudocount_weight)*Pi_true + pseudocount_weight/q*ones(alignment_width, q);

%correct things on the diagonal
scra = eye(q);
for i=1:alignment_width
	for alpha = 1:q
		for beta = 1:q
			Pij(i, i, alpha, beta) = (1.-pseudocount_weight)*Pij_true(i, i, alpha, beta) + pseudocount_weight/q*scra(alpha, beta);
		end
	end
end

end


function C = Compute_C_mat(Pij, Pi, alignment_width, q)
%compute correlation matrix from the matrices of frequencies. 
%Ignore the 1st aa type(=gap) data (it's redundant - freqs sum to one)
%Also reformat to get a square L*(q-1)*L*(q-1) matrix.

C=zeros(alignment_width*(q-1), alignment_width*(q-1));
for i=1:alignment_width
	for j=1:alignment_width
		for alpha=1:q-1
			for beta=1:q-1
				 C(mapkey(i, alpha, q), mapkey(j, beta, q)) = Pij(i, j, alpha+1, beta+1) - Pi(i, alpha+1)*Pi(j, beta+1);
			end
		end
	end
end

end


function W=Compute_DCA_couplings(C, i, j, q)

%get matrix of correlations just for sites i and j
W = zeros(q, q); 

%couplings in 1st line & col, corresponding to gap (residue type 1), are set to ZERO (this
%type of residue is the one for which correlations have been suppressed)
W(2:q, 2:q) = C(mapkey(i, 1:q-1, q), mapkey(j, 1:q-1, q));

end


function W = Change_gauge(W,q)

%perform gauge change to the zero-sum gauge in the block of interest
for i = 1:q
    W(:,i) = W(:,i) - mean(W(:,i));
end
for i = 1:q
    W(i,:) = W(i,:) - mean(W(i,:));
end 

end


function A=mapkey(i, alpha, q)
A = (q-1)*(i-1)+alpha;
end

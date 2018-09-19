function HKRR_energy = Compute_energies(test_seqs,NSeqs,W_store, LengthA, L)
%calculate interaction energy between all pairs of HKs and RRs in test_seqs
%Line: HK; Col: RR.

HKRR_energy = zeros(NSeqs,NSeqs);

for i = 1:NSeqs %to choose the HK
    for j = 1:NSeqs %to choose the RR
        for a = 1:LengthA %sites in HK
            for b = LengthA+1:L %sites in RR
                aa1=test_seqs(i,a); %aa in HK i at site a
                aa2=test_seqs(j,b); %aa in RR j at site b
                HKRR_energy(i,j) = HKRR_energy(i,j) + W_store(a,b,aa1,aa2);
            end
        end
    end
end

end

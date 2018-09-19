function scrambled_alignment = ScrambleSeqs(encoded_focus_alignment, LengthA, table_count_species)
%scrambles the pairings


[N, alignment_width] = size(encoded_focus_alignment);
L=alignment_width-2; % last 2 cols contain species index and initial sequence index!

scrambled_alignment=zeros(N,L+4); %added columns will contain the same data but for the RR


%loop over species
for i=1:size(table_count_species,1)
    
    NSeqs = table_count_species(i,3)-table_count_species(i,2)+1; 
    
    thePerm = randperm(NSeqs); %construct random permutation of numbers from 1 to NSeqs - these will be the indices of the RRs
    
    for j=1:NSeqs
        HKIndex=table_count_species(i,2)+j-1;
        RRIndex=table_count_species(i,2)+thePerm(j)-1;
        scrambled_alignment(HKIndex,1:LengthA)=encoded_focus_alignment(HKIndex,1:LengthA);
        scrambled_alignment(HKIndex,LengthA+1:L)=encoded_focus_alignment(RRIndex,LengthA+1:L);
        scrambled_alignment(HKIndex,L+1:L+2)=encoded_focus_alignment(HKIndex,L+1:L+2);
        scrambled_alignment(HKIndex,L+3:L+4)=encoded_focus_alignment(RRIndex,L+1:L+2);
    end
   
    
end

scrambled_alignment( ~any(scrambled_alignment,2), : ) = [];  %delete rows of zeros (coming from 1st and last seqs, not tabulated in table_count_species)

end

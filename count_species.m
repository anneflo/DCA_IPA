function table_count_species =count_species(encoded_focus_alignment)
%return a table with species id, 1st index in encoded_focus_alignment, last
%index in encoded_focus_alignment

[N, alignment_width] = size(encoded_focus_alignment);
L=alignment_width-2; % last 2 columns contain species index and initial sequence index!

table_count_species=zeros(N,3);
species_id=encoded_focus_alignment(2,L+1); %start at 2nd sequence
iini=2;
count=0;

for i=2:N-1 %run through alignment, avoiding 1st and last sequences (dummies)
    
    species_id_prev=species_id;
    species_id=encoded_focus_alignment(i,L+1);
    
    if species_id~=species_id_prev %a new species is reached for i
        count=count+1;
        ifin=i-1; %last index of previous species
        %save data for previous species
        table_count_species(count,1)=species_id_prev;
        table_count_species(count,2)=iini;
        table_count_species(count,3)=ifin;
        iini=i;
    end
    
end

%save data for last species!
count=count+1;
ifin=N-1;
table_count_species(count,1)=species_id;
table_count_species(count,2)=iini;
table_count_species(count,3)=ifin;

%now suppress the final rows of zeros
table_count_species( ~any(table_count_species,2), : ) = [];

end

function [NUP_alignment, NUP_alignment_headers] = SuppressSpeciesWithOnePair(encoded_focus_alignment, encoded_focus_alignment_headers, table_count_species)
%produce alignment where the species with only one pair are suppressed

N = size(encoded_focus_alignment,1);

ToEl=zeros(N,1);
count=0;

%loop over species
for i=1:size(table_count_species,1)
    
    if table_count_species(i,2)==table_count_species(i,3) %only one pair in this species - needs to be suppressed
        count=count+1;
        ToEl(count,1)=table_count_species(i,2); %final index (i.e. order index in "encoded_focus_alignment") of the sequence to be eliminated       
    end
    
end

%delete rows of zeros at the bottom of ToEl
ToEl( ~any(ToEl,2), : ) = [];

NUP_alignment=encoded_focus_alignment;
NUP_alignment_headers=encoded_focus_alignment_headers;

NUP_alignment(ToEl,:)=[];
NUP_alignment_headers(ToEl)=[];

end

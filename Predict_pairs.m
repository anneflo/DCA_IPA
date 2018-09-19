function Results = Predict_pairs(encoded_focus_alignment, W, LengthA, table_count_species)
%makes pairing predictions on encoded_focus_alignment using DCA scores

[N, alignment_width] = size(encoded_focus_alignment);
L=alignment_width-2; % last 2 cols contain species index and initial sequence index!

%initialize the Results array, used for saving data
Results=zeros(2*N,6);
%col 1: species
%col 2: HK index in initial alignment
%col 3: RR index in initial alignment
%col 4: absolute energy of pairing
%col 5: gap wrt HK
%col 6: gap wrt RR

%total pair counter
totcount = 0; 


%loop over species
for i=1:size(table_count_species,1)
    
    test_seqs = encoded_focus_alignment(table_count_species(i,2):table_count_species(i,3),:);
    NSeqs = table_count_species(i,3)-table_count_species(i,2)+1;
    species_id=table_count_species(i,1);
    
    DeletedInd=zeros(NSeqs,2);
    count=0;

    %now compute the interaction energies of all the HK-RR pairs within the species corresponding to i
    HKRR_energy = Compute_energies(test_seqs,NSeqs,W, LengthA, L);
    HKRR_energy_ini = HKRR_energy;
    
    %now make pairings based on energy within this species (first pair the HK&RR that have the lowest energy etc.)
    
    if NSeqs==1 %trivial case where there is only one pair in this species - make this pair
        
        totcount=totcount+1;
        Results(totcount,1)=species_id;
        Results(totcount,2)=test_seqs(1,L+2); %initial index of sequence
        Results(totcount,3)=Results(totcount,2); %this pairing is right...
        Results(totcount,4)=HKRR_energy;
        Results(totcount,5)=abs(HKRR_energy); %no real gap... consider that absolute energy is gap
        Results(totcount,6)=Results(totcount,5);
        
    else %there is more than one pair in this species - the interesting part starts here!
        
        bigval = max(max(HKRR_energy(:))+10,0);
        val = min(HKRR_energy(:));
        
        
        while val<bigval %not all the HK and RR have been eliminated: there are still pairings to be made, pursue
                        
            
            [r, c] = find(HKRR_energy == val); %indices where the min energy is found => make this pair          
            
            %compute energy gaps associated to these possible pairings
            theGaps=zeros(size(r,1),5);
            theGaps(:,1)=r;
            theGaps(:,2)=c;           
            
            for myind=1:size(r,1)
                
                %first compute gap wrt HK, i.e. within row
                theLine=HKRR_energy(r(myind),:);
                theLine(theLine == val)=bigval;
                val2=min(theLine);%next smallest energy of a coupling for this HK (behind val)
                if val2==bigval %there is only one RR left that can pair with this HK. Gap notion is ill-defined.
                    theGaps(myind,3)=-1;
                else
                    theGaps(myind,3)=val2-val;
                    better=HKRR_energy_ini(r(myind),DeletedInd(1:count,2));
                    better=better(better<val);
                    if ~isempty(better) %the best match(es) have been deleted by previous matches - penalize
                        theGaps(myind,3)=(val2-val)./(length(better)+1);
                    end
                end
                
                %then compute gap wrt RR, i.e. within col
                theCol=HKRR_energy(:,c(myind));
                theCol(theCol == val)=bigval;
                val2=min(theCol);%next smallest energy of a coupling for this RR (behind val)
                if val2==bigval %there is only one HK left that can pair with this RR. Gap notion is ill-defined.
                    theGaps(myind,4)=-1;
                else
                    theGaps(myind,4)=val2-val;
                    better=HKRR_energy_ini(DeletedInd(1:count,1),c(myind));
                    better=better(better<val);
                    if ~isempty(better) %the best match(es) have been deleted by previous matches - penalize
                        theGaps(myind,4)=(val2-val)./(length(better)+1);
                    end
                end
                
                %also store min gap
                theGaps(:,5)=min(theGaps(:,3),theGaps(:,4));
                
            end
            
            
            if size(r,1)>1
                %rank the ex-aequo best pairs by min gap
                auxGap=find(theGaps(:,5)==min(theGaps(:,5)));
                if size(auxGap,1)==1
                    theGaps=theGaps(auxGap,:);
                else %pick a random pair
                    theGaps=theGaps(randi(size(r,1)),:);
                end
            end
                        
            %make one pair and then move on                                                     
            totcount=totcount+1;                                                  
            count=count+1;

            Results(totcount,1)=species_id;
            Results(totcount,2)=test_seqs(theGaps(1,1),L+2); %initial index of HK sequence (HK: line)
            Results(totcount,3)=test_seqs(theGaps(1,2),L+2); %initial index of RR sequence (RR: col)
            Results(totcount,4)=val; %absolute energy of the pairing - was the min of the HKRR_energy matrix
            Results(totcount,5)=theGaps(1,3); %gap within row - wrt this HK
            Results(totcount,6)=theGaps(1,4); %gap within col - wrt this RR  

            DeletedInd(count,1)=theGaps(1,1);
            DeletedInd(count,2)=theGaps(1,2);

            %since we only want one-to-one pairings, the HK and RR can't interact with anything: suppress them from further consideration
            HKRR_energy(theGaps(1,1),:) = bigval*ones(1,NSeqs); %replace the row by bigval to effectively suppress HK from further consideration
            HKRR_energy(:,theGaps(1,2)) = bigval*ones(NSeqs,1); %replace the col by bigval to effectively suppress RR from further consideration               
            
            %now update val for next round!
            val = min(HKRR_energy(:));
            
        end
        
        %Almost done with this species - all pairings have been made
        %Just fix the undefined gaps
        aux=Results(Results(:,1)==species_id,5);
        aux(aux==-1)=Inf;
        if min(aux)~=Inf 
            Results(Results(:,1)==species_id & Results(:,5)==-1,5)=min(aux);
        else
            Results(Results(:,1)==species_id & Results(:,5)==-1,5)=0;
        end
        aux=Results(Results(:,1)==species_id,6);
        aux(aux==-1)=Inf;
        if min(aux)~=Inf 
            Results(Results(:,1)==species_id & Results(:,6)==-1,6)=min(aux);
        else
            Results(Results(:,1)==species_id & Results(:,6)==-1,6)=0;
        end
        
        
    end
    
end

%delete rows of zeros
Results( ~any(Results,2), : ) = [];

end

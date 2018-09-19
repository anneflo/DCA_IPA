%This is the main code to run the DCA-IPA on the standard HK-RR dataset.
clear all
close all hidden

%set parameters
Nincrement = 1600;
LengthA = 64;%length of first protein (here the HK)

%read data files
msa_fasta_filename = 'Standard_HKRR_dataset.fasta'; %sequence data file
load SpeciesNumbering_Standard_HKRR_dataset; %read in SpeciesNumbering_extr

%read sequences, adding species number in L+1 and sequence number in L+2
%L is the full length of concatenated sequences, without supplementary indicators such as species and initial index
[encoded_focus_alignment, encoded_focus_alignment_headers, L] = readAlignment_and_NumberSpecies(msa_fasta_filename,SpeciesNumbering_extr);
%L is alignment_width (full length of seqs, without supplementary indicators such as species and initial index)
disp(L)

%suppress species with one pair
table_count_species =count_species(encoded_focus_alignment);
[encoded_focus_alignment, encoded_focus_alignment_headers] = SuppressSpeciesWithOnePair(encoded_focus_alignment, encoded_focus_alignment_headers, table_count_species);
N = size(encoded_focus_alignment,1); %number of sequences
disp(N)
%tabulate species and sequences within species
table_count_species =count_species(encoded_focus_alignment);

%number of rounds (last one -> all sequences are in the training set)
Nrounds=ceil(N./Nincrement+1); 
disp(Nrounds)

%start from random within-species pairings: scramble the pairings for this.
encoded_training_alignment = ScrambleSeqs(encoded_focus_alignment, LengthA, table_count_species);
%save the species and initial indices of the sequences in the scrambled alignment we start from
filename=strcat('Res/IniScrambling_Ninc',num2str(Nincrement),'.txt');
dlmwrite(filename,encoded_training_alignment(:,L+1:L+4),'delimiter','\t')
%in the training set, discard extra indices (species index, initial sequence index)
encoded_training_alignment(:,L+1:L+4)=[];

%initialize
NSeqs_new=0;
Output=zeros(Nrounds,6); 
%Output matrix:
%Each row corresponds to an iteration of the DCA-IPA.
%col 1: number of sequences NSeqs in concatenated alignment used as training set
%col 2: effective number of sequences Meff in concatenated alignment used as training set
%col 3: number of TP pairs
%col 4: number of FP pairs
%col 5: number of TP pairs in concatenated alignment used as training set
%col 6: number of FP pairs in concatenated alignment used as training set


%%

for rounds=1:Nrounds %iterate the process until all sequences end up in the training set

    disp(rounds)
    
    if rounds>1 
        
        %update the training set by adding in the pairings with largest energy gaps made at previous round

        Results(:,6)=min(Results(:,5),Results(:,6)); %compute min of 2 gap scores
        Results=sortrows(Results,-6); %sort by gap score
        
        %number of sequences that will be added to form the training set for this round
        NSeqs_new = NSeqs_new + Nincrement; 
        if NSeqs_new>=size(Results,1)
            NSeqs_new=size(Results,1); %for the last round, all paired sequences will be in the training set
        end
        
        %save to Output the number of TP or FP in the training set
        Output(rounds,5)=size(Results(Results(1:NSeqs_new,2)==Results(1:NSeqs_new,3)),1);
        Output(rounds,6)=size(Results(Results(1:NSeqs_new,2)~=Results(1:NSeqs_new,3)),1);

        %construct new training set
        newseqs=zeros(NSeqs_new,L);
        for i=1:NSeqs_new
            newseqs(i,1:LengthA) = encoded_focus_alignment(encoded_focus_alignment(:,L+2)==Results(i,2), 1:LengthA );
            newseqs(i,LengthA+1:L) = encoded_focus_alignment(encoded_focus_alignment(:,L+2)==Results(i,3), LengthA+1:L );  
        end
        encoded_training_alignment = newseqs; 
        
    end
    
    %construct model from training set
    [W, Meff] = Compute_DCA_model(encoded_training_alignment);
 
    %compute pairings and gap scores for all pairs
    Results = Predict_pairs(encoded_focus_alignment, W, LengthA, table_count_species);  
    
    %save the results
    Output(rounds,1)=NSeqs_new;
    Output(rounds,2)=Meff;
    Output(rounds,3)=size(Results(Results(:,2)==Results(:,3)),1);
    Output(rounds,4)=size(Results(Results(:,2)~=Results(:,3)),1);
    
end


%%
%save Output matrix
filename=strcat('Res/TP_data_Ninc',num2str(Nincrement),'.txt');
dlmwrite(filename,Output,'delimiter','\t')

%save the final pairs made and their scores
filename=strcat('Res/Resf_Ninc',num2str(Nincrement),'.txt');
dlmwrite(filename,Results,'delimiter','\t')

function [encoded_focus_alignment, encoded_focus_alignment_headers, alignment_width] = readAlignment_and_NumberSpecies(msa_fasta_filename,SpeciesNumbering)

%THIS ASSUMES THAT THE FIRST SEQUENCE IS A REFERENCE SEQUENCE
%and that the last one is a dummy

% read in alignment. 
disp(msa_fasta_filename)
full_alignment = fastaread(msa_fasta_filename);
alignment_height = size(full_alignment, 1);
alignment_width = size(full_alignment(1).Sequence, 2);

disp(alignment_height)

%initialize
encoded_focus_alignment = zeros(alignment_height,alignment_width+2); %2 extra cols for: species ID; line in initial alignment
encoded_focus_alignment_headers=cell(alignment_height,1);

%letter to number mapping (defined below) to resolve ambiguities etc.
letter2number_map = create_letter2number_map();


for index=1:alignment_height %loop over sequences in the alignment 
    
	encoded_focus_alignment_row = letter2number_map(full_alignment(index).Sequence); %convert sequence to line vector of numbers    
	encoded_focus_alignment(index,1:alignment_width) = encoded_focus_alignment_row; %append sequence to "encoded_focus_alignment"
    encoded_focus_alignment_headers{index} = full_alignment(index).Header; 
    
    if index>1 && index<alignment_height %ref sequence doesn't have the usual type of header; last may be from different species in this extracted file...
        %get species name from the sequence header
        uniprot_range_line = full_alignment(index).Header;
        Sepposition = strfind(uniprot_range_line,'|');
        species_id = uniprot_range_line(Sepposition(1)+1:Sepposition(2)-1);  
        %use SpeciesNumberingTable to convert this name to a number
        encoded_focus_alignment(index,alignment_width+1) = SpeciesNumbering{strcmp(SpeciesNumbering(:,2),species_id),1};
    end
    
    if index<alignment_height
        encoded_focus_alignment(index,alignment_width+2) = index;
    else
        encoded_focus_alignment(index,alignment_width+1) = NaN;
        encoded_focus_alignment(index,alignment_width+2) = NaN;
    end
        
end

%now suppress any final rows that are empty (=contain only zeros) in the alignment
encoded_focus_alignment( ~any(encoded_focus_alignment,2), : ) = [];

end



%%auxiliary function definitions

function letter2number_map = create_letter2number_map()
%characters are indexed by decimal ASCII codes (which range from 0 to 255, with - ->45 and A ->65 etc.)
letter2number_map(256) = 0; %initialize all bytes to 0 (1 line and 256 cols of 0s)
letter2number_map('-') = 1;
letter2number_map('A') = 2;
letter2number_map('C') = 3;
letter2number_map('D') = 4;
letter2number_map('E') = 5;
letter2number_map('F') = 6;
letter2number_map('G') = 7;
letter2number_map('H') = 8;
letter2number_map('I') = 9;
letter2number_map('K') = 10;
letter2number_map('L') = 11;
letter2number_map('M') = 12;
letter2number_map('N') = 13;
letter2number_map('P') = 14;
letter2number_map('Q') = 15;
letter2number_map('R') = 16;
letter2number_map('S') = 17;
letter2number_map('T') = 18;
letter2number_map('V') = 19;
letter2number_map('W') = 20;
letter2number_map('Y') = 21;
letter2number_map('B') = -1; %ambiguous : skip sequences containing these
letter2number_map('Z') = -1; %ambiguous : skip sequences containing these
letter2number_map('J') = -1; %ambiguous : skip sequences containing these
letter2number_map('X') = -1; %ambiguous : skip sequences containing these
letter2number_map('U') = -1; %non-standard : skip sequences containing these
letter2number_map('O') = -1; %non-standard : skip sequences containing these
letter2number_map('a') = -2; %non-conserved: skip in seq of interest
letter2number_map('c') = -2; %non-conserved: skip in seq of interest
letter2number_map('d') = -2; %non-conserved: skip in seq of interest
letter2number_map('e') = -2; %non-conserved: skip in seq of interest
letter2number_map('f') = -2; %non-conserved: skip in seq of interest
letter2number_map('g') = -2; %non-conserved: skip in seq of interest
letter2number_map('h') = -2; %non-conserved: skip in seq of interest
letter2number_map('i') = -2; %non-conserved: skip in seq of interest
letter2number_map('k') = -2; %non-conserved: skip in seq of interest
letter2number_map('l') = -2; %non-conserved: skip in seq of interest
letter2number_map('m') = -2; %non-conserved: skip in seq of interest
letter2number_map('n') = -2; %non-conserved: skip in seq of interest
letter2number_map('p') = -2; %non-conserved: skip in seq of interest
letter2number_map('q') = -2; %non-conserved: skip in seq of interest
letter2number_map('r') = -2; %non-conserved: skip in seq of interest
letter2number_map('s') = -2; %non-conserved: skip in seq of interest
letter2number_map('t') = -2; %non-conserved: skip in seq of interest
letter2number_map('v') = -2; %non-conserved: skip in seq of interest
letter2number_map('w') = -2; %non-conserved: skip in seq of interest
letter2number_map('y') = -2; %non-conserved: skip in seq of interest
letter2number_map('b') = -2; %non-conserved: skip in seq of interest
letter2number_map('z') = -2; %non-conserved: skip in seq of interest
letter2number_map('j') = -2; %non-conserved: skip in seq of interest
letter2number_map('x') = -2; %non-conserved: skip in seq of interest
letter2number_map('u') = -2; %non-conserved: skip in seq of interest
letter2number_map('o') = -2; %non-conserved: skip in seq of interest
letter2number_map('.') = -3; %non-conserved: skip in seq of interest, do not advance position
end


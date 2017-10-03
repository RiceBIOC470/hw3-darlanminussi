% Darlan Conterno Minussi
%HW3

%% Problem 1 - Smith-Waterman alignment
% Consider two sequences 'GTAATCC' and 'GTATCCG'

seq1 = 'GTAATCC';
seq2 = 'GTATCCG';

matchval = 2;
mismatchval = -1;
gapval = -1;

ofdiag = ones(6)-eye(6);
S = matchval*eye(6)+mismatchval*ofdiag;

[score1, align1, start1] = swalign(seq1, seq2,'Alphabet','nt', 'ScoringMatrix', S, 'Showscore', true);

% Construct the scoring matrix for this with the parameters:
% match value = 2, mismatch value = -1, and gap penalty = -1. Use your
% solution to get the optimal alignment. If you prefer, it is acceptable to do this with
% pencil and paper, you can then take a snapshot of your solution and
% include it in your repository. 

%% Problem 2 - using the NCBI databases and sequence alignments

% Erk proteins are critical signal transducers of MAP kinase signaling.
% Accessions numbers for ERK1 (also called MAPK3) and ERK2 (also called MAPK1) human mRNA are NM_002746 and
% NM_002745, respectively. 

erk1 = getgenbank('NM_002746', 'SequenceOnly', true);
erk2 = getgenbank('NM_002745', 'SequenceOnly', true);

% Part 1. Perform an alignment of the coding DNA sequences of ERK1 and
% ERK2. What fraction of base pairs in ERK1 can align to ERK2? 

[score2, align2, start2] = swalign(erk1, erk2, 'Alphabet','nt', 'Showscore', true);

count_align2 = countAligned(align2);
disp(count_align2);


% Part2. Perform an alignment of the aminoacid sequences of ERK1 and ERK2.
% What fraction of amino acids align?

erk1pt = getgenbank('NM_002746');
erk2pt = getgenbank('NM_002745');

[score2, align2, start2] = swalign(erk1pt.CDS.translation, erk2pt.CDS.translation,'Alphabet', 'aa');

countAlignedaa(align2);


% Part 3.  Use the NCBI tools to get mRNA sequences for the mouse genes ERK1 and
% ERK2 and align both the coding DNA sequences and protein sequences to the
% human versions. How similar are they? 

erk1mouse = getgenbank('NM_011952');
erk2mouse = getgenbank('NM_011949');


% Aligning mRNA erk1 mouse with human erk1 mRNA
[score3, align3, start3] = swalign(erk1, erk1mouse.Sequence, 'Alphabet', 'nt');
count_align3 = countAligned(align3);
disp(count_align3);

% Aligning hERK1 protein with mERK1 protein
[score4, align4, start4] = swalign(erk1pt.CDS.translation, erk1mouse.CDS.translation, 'Alphabet', 'aa');
countAlignedaa(align4);

% Aligning erk2 mouse mRNA with human erk2 mRNA
[score5, align5, start5] = swalign(erk2, erk2mouse.Sequence, 'Alphabet', 'nt');
count_align5 = countAligned(align5);
disp(count_align5);

% Aligning hERK2 protein with mERK2 protein
[score6, align6, start6] = swalign(erk2pt.CDS.translation, erk2mouse.CDS.translation, 'Alphabet', 'aa');
countAlignedaa(align6);

%% Problem 3: using blast tools programatically

% Part 1. Write a function that takes an NCBI accession number and a number N as input and
% returns a cell array of the accession numbers for the top N blast hits. 

[all, no_h, only_h] = blastncbi_hits('NM_011952', 6);
disp(all);

% Part 2. Write a function that takes an accession number as input, calls your function 
% from part 1, and returns two outputs - the closest scoring match in human DNA/RNA and the
% closest non-human match. Hint: see the "Source" or "SourceOrganism" field in the data
% returned by getgenbank. Make sure your function does something sensible
% if nothing human is found. 

[all, no_h, only_h] = blastncbi_hits('NM_011952', 1);
disp(all);
disp(no_h);

% Part 3. Choose at least one gene from the human genome and one gene from
% another organism and run your code from part 2 on the relevant accession
% numbers. Comment on the results. 

% human TNFAIP3
[all_tnfaip3, no_h_tnfaip3, only_h_tnfaip3] = blastncbi_hits('NM_001270508', 6);
disp(all_tnfaip3);
disp(no_h_tnfaip3);
disp(only_h_tnfaip3);


% Mus musculus TNFAIP3

[all_tnfaip3mus, no_h_tnfaip3mus, only_h_tnfaip3mus] = blastncbi_hits('NM_009397', 6);
disp(all_tnfaip3mus);
disp(no_h_tnfaip3mus);
disp(only_h_tnfaip3mus);

% Here we compared results from a search of the same gene TNFAIP3 which
% encodes for a ubiquitin editing protein.
% When we search for the human mRNA as expected, the first results are
% transcript variants present in humans, and the other closest results are
% from commom chimpanzee whereas in the search for Mus musculus the top
% scoring results are from transcript variants observed in mice and ESTs in
% bone marrow macrophage of mice
% interestingly, when we remove all the scoring results from the species
% Homo sapiens, TNAFAIP3 from Mus musculus does not appear in the top 6
% observed results
% As denoted by the message, no human sequences were shown in the sequences
% returned from blast when searching for the mouse TNFAIP3 gene.



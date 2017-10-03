function [accession_hits, no_human, only_human] = blastncbi_hits(accession_number, N)

% INPUT: Accession number of the desired transcript or protein to be
% searched in blast
% N number of results to return
% OUTPUT: cell array of the accession numbers for the top N blast hits.

genbank_info = getgenbank(accession_number);

[RID, RT0E] = blastncbi(genbank_info.Sequence, 'blastn');

blastdata = getblast(RID, 'WaitTime', RT0E + 10);

blast_table = struct2table(blastdata.Hits);
blast_cell = table2cell(blast_table);

% removing Homo Sapiens
blast_table_noh = blast_table(cellfun(@isempty, regexp(blast_table.Name, 'Homo sapiens')), :);

% keeping only Homo Sapiens
blast_table_only_h = blast_table(~cellfun(@isempty, regexp(blast_table.Name, 'Homo sapiens')), :);

accession_hits = {};
no_human = {};
only_human = {};

for i = 1:N
    accession_hits{i,1} = string(blast_table.Name(i));
end


if isempty(blast_table_only_h)
    disp('No results were from Homo Sapiens');
else
    if N > size(blast_table_only_h,1)
        P = size(blast_table_only_h,1);
        for i = 1:P
            only_human{i,1} = string(blast_table_only_h.Name(i));
        end
    end
end

for i = 1:N
    no_human{i,1} = string(blast_table_noh.Name(i));
end




end


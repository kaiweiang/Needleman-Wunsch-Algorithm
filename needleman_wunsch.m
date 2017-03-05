function [score, alignment1, aligner, alignment2] = needleman_wunsch(seq1, seq2, match, mismatch, gap)
% needleman_wunsch Perform the Needleman-Wunsch algorithm.
%
% Input variables:
% seq1: string containing the first sequence
% seq2: string containing the second sequence
% match: numeric value which represents the score for match
% mismatch: numeric value which represents the penalty for mismatch
% gap: numeric value which represents the penalty for gap
%
% Output variables:
% score: numeric value containing the alignment score
% alignment1: string containing the alignment (including gaps) of the first
% sequence
% aligner: string containing symbols (pipe |) which aligns the 2 sequences
% alignment2: string containing the alignment (including gaps) of the
% second sequence
%
% Output example:
% AV--TNAGQLV---S  <--- alignment1
% ||  || |||    |  <--- aligner
% AQVSTN-GQL-AQVT  <--- alignment2
%
%
%%%%%%%%%%%%%% YOUR CODE STARTS HERE
seq1_len = length(seq1);
seq2_len = length(seq2);
F = zeros(seq2_len+1, seq1_len+1); %preallocate memory for the scoring matrix
F(1,2:end) = gap*(1:seq1_len); 
F(2:end,1) = gap*(1:seq2_len); 
TB = repmat('x',size(F)); %preallocate memory for the trace back matrix
TB(1:end,1) = repmat('u',size(TB,1),1); 
TB(1,1:end) = repmat('l',1,size(TB,2));  
TB(1,1) = '0';
keys = {1,2,3};
values = {'d','u','l'}; %d= diagonal; u = up; l = left
dict = containers.Map(keys,values);

%Needleman-Wunsch Algorithm
for i = 2:size(F,1) %row
    for j = 2:size(F,2) %col
        
        %case 1
        if seq2(i-1) == seq1(j-1)
            case1 = F(i-1,j-1) + match;
        else
            case1 = F(i-1,j-1) + mismatch;
        end
        
        %case 2
        case2 = F(i-1,j) + gap;
        
        %case 3
        case3 = F(i,j-1) + gap;
        
        [max_score,max_score_idx] = max([case1 case2 case3]);
        F(i,j) = max_score;
        TB(i,j) = dict(max_score_idx);
        
    end
end

%Trace back
i = size(TB,1); j = size(TB,2);
tb=[]; %tb contains the optimal path found from the trace back matrix
while i>1 || j>1
    tb = [tb TB(i,j)];
    switch TB(i,j)
        case 'd'
            i = i-1; j = j-1;
        case 'u'
            i = i-1;
        case 'l'
            j = j-1;
    end
end
tb = fliplr(tb); 

%Add gaps to the sequences
seq1_gap_idx = find(tb=='u');
logical_idx = false(1,length(seq1)+length(seq1_gap_idx));
logical_idx(seq1_gap_idx) = true;
new_seq1 = nan(size(logical_idx));
new_seq1(~logical_idx) = seq1;
new_seq1 = char(new_seq1);
new_seq1(logical_idx) = '-';

seq2_gap_idx = find(tb=='l');
logical_idx = false(1,length(seq2)+length(seq2_gap_idx));
logical_idx(seq2_gap_idx) = true;
new_seq2 = nan(size(logical_idx));
new_seq2(~logical_idx) = seq2;
new_seq2 = char(new_seq2);
new_seq2(logical_idx) = '-';

aligner = repmat(' ',size(new_seq1));
aligner(new_seq1 ~= '-' & new_seq2 ~= '-') = '|'; % add the vertical line for each match and mismatch
alignment1 = new_seq1;
alignment2 = new_seq2;
score = F(end,end);

end
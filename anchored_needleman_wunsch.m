function [score, alignment1, aligner, alignment2] = anchored_needleman_wunsch(seq1, seq2, match, mismatch, gap, matched_regions)
% anchored_needleman_wunsch Perform the anchored version of the 
% Needleman-Wunsch algorithm.
%
% Input variables:
% seq1: string containing the first sequence
% seq2: string containing the second sequence
% match: numeric value which represents the score for match
% mismatch: numeric value which represents the penalty for mismatch
% gap: numeric value which represents the penalty for gap
% matched_regions: contains regions which are know to be aligned. The
% regions may contain matches and/or mismatches. The first column contains
% the starting position of the first sequence; the second column contains
% the ending position of the first sequence; the third column contains the
% starting position of the second sequence; and the forth column contains
% the ending position of the second sequence. E.g.: a row containing [17,
% 20, 76, 79] means that seq1(17:20) should be aligned to seq2(76:79),
% containing only matches and/or mismatches in this region.
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
num_matched_regions = size(matched_regions,1);
SCORE = zeros(1,num_matched_regions*2+1); %preallocate memory for the regional score vector 

%get the score of the first unaligned region
[sc,a1,alig,a2] = needleman_wunsch(seq1(1:matched_regions(1,1)-1), seq2(1:matched_regions(1,3)-1), match, mismatch, gap);
SCORE(1) = sc;
alignment1 = a1; aligner = alig; alignment2 = a2;

count = 2; %index for the SCORE
for n = 1:num_matched_regions
    
    %matched region
    aligned_seq1 = seq1(matched_regions(n,1):matched_regions(n,2));
    aligned_seq2 = seq2(matched_regions(n,3):matched_regions(n,4));
    sc = sum(aligned_seq1 == aligned_seq2)*match;
    sc = sc + (sum(aligned_seq1 ~= aligned_seq2)*mismatch);
    alig = repmat('|',size(aligned_seq1));
    SCORE(count) = sc;
    alignment1 = strcat(alignment1, aligned_seq1);
    aligner = [aligner, alig];
    alignment2 = strcat(alignment2, aligned_seq2);
    
    %non matched region
    if n ~= num_matched_regions       
        unaligned_seq1 = seq1(matched_regions(n,2)+1:matched_regions(n+1,1)-1);
        unaligned_seq2 = seq2(matched_regions(n,4)+1:matched_regions(n+1,3)-1);
    else
        unaligned_seq1 = seq1(matched_regions(n,2)+1:end);
        unaligned_seq2 = seq2(matched_regions(n,4)+1:end);
    end
    
    [sc,a1,alig,a2] = needleman_wunsch(unaligned_seq1, unaligned_seq2, match, mismatch, gap);
    alignment1 = strcat(alignment1, a1);
    aligner = [aligner, alig];
    alignment2 = strcat(alignment2, a2);
    SCORE(count+1) = sc;
    count = count + 2;    
end
score = sum(SCORE); %total score

end


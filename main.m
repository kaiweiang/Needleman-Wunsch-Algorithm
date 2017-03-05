%% load sequences
seq1 = fastaread('sequences/Human_HOX.fa'); 
seq1 = seq1.Sequence;
seq2 = fastaread('sequences/Fly_HOX.fa');
seq2 = seq2.Sequence;

%% penalties
gap = -2;
match = 1;
mismatch = -3;

%% anchored needleman wunsch
matched_regions = dlmread('sequences/Match_HOX.txt');
[score, alignment1, aligner, alignment2] = anchored_needleman_wunsch(seq1, seq2, match, mismatch, gap, matched_regions);

%% perform alignment of random sequences 10,000 times
acc_scores = zeros(1,1e4);
seq1_len = length(seq1);
seq2_len = length(seq2);
parfor i=1:100%1e4
    % generate random sequence 1 and 2 using seq1 and seq2
    seq1_rand = seq1(randperm(seq1_len));
    seq2_rand = seq2(randperm(seq2_len));

    % call anchored_needleman_wunsch and store scores into acc_scores
    [scr,~,~,~] = anchored_needleman_wunsch(seq1_rand, seq2_rand, match, mismatch, gap, matched_regions);    
    acc_scores(i) = scr;
end

%% plot random scores and actual score (uncomment below to plot)
figure1 = figure;
axes1 = axes('Parent',figure1);
histogram(acc_scores)
hold on
histogram(score)
grid on; grid minor;
legend({'Random Sequences', 'Alignment Score'});
xlabel('Score'); ylabel('Frequency');
set(axes1,'FontSize',14);

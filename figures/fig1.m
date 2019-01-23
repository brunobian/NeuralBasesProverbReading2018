%% Load
clear all
close all
load('../data/exp')
addpath('../functions/')

%% Fig 1 (b)(c) - Methods  - Corpus
close all

catgram = [exp.catgram];
disp(['Total words: ' num2str(length(catgram))])
disp(['Total content words: ' num2str(sum(ismember(catgram, 'anv')))])

words = {};
catgram = [];
allWords = [];

for i=1:120
    sntc = strsplit(exp(i).sentence, ' ');
    catgramSntc = exp(i).catgram;
    for j=1:length(sntc)
        
        pal = sntc{j};
        cat = catgramSntc(j);
        tmp = strfind(words, pal);
        
        allWords = [allWords, pal];
        if ~any(cellfun(@(x) any(x), tmp))
            words = [words, pal];
            catgram = [catgram, cat];
        end
    end
end
disp(['Unique words: ' num2str(length(words))])
disp(['Unique content words: ' num2str(sum(ismember(catgram, 'anv')))])

pred = [exp.logitPred];
disp(['mean pred proverbss: ' num2str(mean(pred([exp.typeWrd] == 0)))])
disp(['SD pred   proverbss: ' num2str(std(pred([exp.typeWrd] == 0)))])
disp(['mean pred No proverbss: ' num2str(mean(pred([exp.typeWrd] ~= 0)))])
disp(['SD pred   No proverbss: ' num2str(std(pred([exp.typeWrd] ~= 0)))])


% Logicals
provs  = [exp.typeWrd] == 0 & [exp.filter];
nProvs = [exp.typeWrd] ~= 0 & [exp.filter];

posRelRP = [exp.posRelRP];
pred = [exp.logitPred];

% Dist de Pred
[yProv, xProv]   = hist(pred(provs));  yProv = yProv / sum(provs);
[ynProv, xnProv] = hist(pred(nProvs)); ynProv = ynProv / sum(nProvs);
figure();clf
set(gcf,'Color','w')
hold on
    plot(xProv, yProv, 'r', 'LineWidth', 2)    
    plot(xnProv, ynProv, 'LineWidth', 2)
    xlabel('logit(Pred)')
    ylabel('Proportion of words')
    legend({'Memory-enconded','Common Sentences'})

    
% Mean MJ
unPalNumRelMj = unique(posRelRP(provs));
predRel =[];
for i = 1:length(unPalNumRelMj)
    thisPals = posRelRP == unPalNumRelMj(i);
    predRel = [predRel median(pred(thisPals))];
end

figure();clf
set(gcf,'Color','w')
hold on
    scatter(posRelRP(provs), pred(provs))
    plot(unPalNumRelMj, predRel, 'r', 'LineWidth', 2)
    xlim([-4 7])
    xlabel('Position relative to MJ')
    ylabel('logit(Pred)')
    legend({'raw','median'}, 'location', 'SouthEast')

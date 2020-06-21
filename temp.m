iterations = 10000;
trials = 90;
pa = .1;
pb = .1;
minInterval = 3;


[empA, empB] = deal(nan(1,iterations));

for j = 1:iterations
    trialType = nan(1,trials);
    count = 0;

    % adjust probabilities
    paAdj = - pa / (pa*minInterval + pb*minInterval - 1);
    pbAdj = - pb / (pa*minInterval + pb*minInterval - 1);

    for i = 1:trials
        r = rand(1);
        if r<paAdj && count>=minInterval
            trialType(i) = 1;
            count = 0;
        elseif r<(paAdj+pbAdj)  && count>=minInterval
            trialType(i) = 2;
            count = 0;
        else
            trialType(i) = 3;
            count = count + 1;
        end
    end

%     fprintf('\ntrial proportions...\n')
%     fprintf('a: %.3f\n', mean(trialType==1))
%     fprintf('b: %.3f\n', mean(trialType==2))
%     fprintf('c: %.3f\n', mean(trialType==3))
    empA(j) = mean(trialType==1);
    empB(j) = mean(trialType==2);
end
disp('all done!')

%%
figure; hold on;
histogram(empA, 'Normalization', 'probability')
histogram(empB, 'Normalization', 'probability')






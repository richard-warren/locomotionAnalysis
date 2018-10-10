


sessions = {'180917_002', '180920_002', '180922_001', '181001_002', '181002_002', '181003_002', '181004_003'};

for i = 1:length(sessions)
    try
        plotRecordingSummary(sessions{i})
    catch
        fprintf('%s: PROBLEM ANALYZING SESSION!\n', sessions{i})
    end
end
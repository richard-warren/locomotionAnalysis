function sessionDvs = getSessionDvs(speedAvoidanceData, kinData, dvs)

% given speedAvoidanceData and kinData structures (which contain trial to
% trial information for multiple sessions), computes dependent measures
% averaged across each session // only computes dependent measures listed
% in dvs cell array

keyboard

% initializations
kinSessions = unique({kinData.session});
speedAvoidanceSessions = unique({speedAvoidanceData.session});
if isequal(kinSessions, speedAvoidanceSessions); sessions = kinSessions; else; disp('WARNING! different sessions detected in kinData and speedAvoidanceData!'); end


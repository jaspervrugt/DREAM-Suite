function prt_progress(DREAMPar,N)
% Print progress

persistent count ct
if isempty(count)
    count = 0; ct = 0; fprintf('\n');
else
    ct = ct + DREAMPar.CPU; 
    % It is possible that ct does not get to N because we have a remaining
    % balance unequal to DREAMPar.CPU - fix this for last iteration
end
fprintf(1, repmat('\b',1,count)); % delete line before
count = fprintf('Posterior simulation, %% done: %3.2f',100 * ( ct/N ) );
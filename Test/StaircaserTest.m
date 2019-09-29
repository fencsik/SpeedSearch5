function StaircaserTest

rand('twister', 100 * sum(clock));
thresholdMean = 80;
thresholdSD = 10;

id = Staircaser('Create', 1, 10, 50, [-5 10 0], 2, 3, [0 100]);
[progress, stepsize] = Staircaser('Progress', id);
fprintf('\nStaircase progress stepsize = %0.3f\n\n', stepsize);

done = 0;
trial = 0;

while ~done
    trial = trial + 1;
    a = thresholdMean + thresholdSD .* randn;
    [success, val, track] = Staircaser('StartTrial', id);
    if mod(trial, 12) == 0
        resp = 3;
        respString = 'MAYBE';
    elseif mod(trial, 17) == 0
        resp = 0;
        respString = 'BAD';
    elseif val >= a
        resp = 1;
        respString = 'CORRECT';
    else
        resp = 2;
        respString = 'ERROR';
    end
    [success, done] = Staircaser('EndTrial', id, resp);
    progress = Staircaser('Progress', id);
    fprintf('TRIAL %d: VALUE %d, %s (progress = %0.3f)\n\n', trial, val, respString, progress);
end

val = Staircaser('FinalValue', id);
rev = Staircaser('GetReversals', id)

fprintf('Final Staircase Value = %0.3f\n\n', val);

% Staircaser('Plot', id);

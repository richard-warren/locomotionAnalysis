function times = flirTimeStampDecoderFLIR(timeStamps)

    % converts metadata timestamps recorded in bonsai from point grey camera into seconds
    %
    % input:     timeStamps: metadata timestamps recorded in bonsai
    % output:    times: seconds, normalized such that first time is 0
    
        
    % extract first cycle (first 7 bits)
    cycle1 = bitshift(timeStamps, -25);

    % extract second cycle (following 13 bits)
    cycle2 = bitand(bitshift(timeStamps, -12), 8191) / 8000; % 8191 masks all but but last 13 bits
    
    % account for overflows in counter
    times = cycle1 + cycle2;
    overflows = [0; diff(times) < 0];
    times = times + (cumsum(overflows) * 128);
    
    % offset such that first time is zero
    times = times - min(times);
    
end
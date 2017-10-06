function times = timestampDecoder(timeStamps)

    % converts metadata timestamps recorded in bonsai from point grey camera into seconds
    %
    % input:     timeStamps: metadata timestamps recorded in bonsai
    % output:    times: seconds, normalized such that first time is 0
    
    
    
    % extract first cycle (following 13 bits)
    cycle1 = de2bi(bitshift(timeStamps, -25), 'left-msb');
    cycle1 = bi2de(cycle1, 'left-msb');

    % extract first cycle (first 7 bits)
    cycle2 = bitand( de2bi(bitshift(timeStamps, -12), 'left-msb'),  [zeros(1,7) ones(1,13)] );
    cycle2 = bi2de(cycle2, 'left-msb') / 8000;
    
    % account for overflows in counter
    times = cycle1 + cycle2;
    overflows = [0; diff(times) < 0];
    times = times + (cumsum(overflows) * 128);
    
    % offset such that first time is zero
    times = times - min(times);
    
end
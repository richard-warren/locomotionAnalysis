


smps = 1000;
bits = 8;

values = round((sin(linspace(3*pi/2, 3*pi/2+2*pi, smps))/2+.5) * (2^bits-1));
close all; figure; plot(values);
clipboard('copy', sprintf('%i,', values)) % vopy values to clipboard

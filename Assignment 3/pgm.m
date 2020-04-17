function [pgm_out] = pgmest(x)
    N = length(x);
    pgm_out = (1/N)*abs(fft(x)).^2; 
end
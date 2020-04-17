function [p_out] = pgm(x)
    N = length(x);
    p_out = (1/N)*abs(fft(x)).^2; 
end
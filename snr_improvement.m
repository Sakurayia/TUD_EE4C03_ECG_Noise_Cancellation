function y = snr_improvement(original, mixed, reconstructed)

snr_original = snr(original, mixed - original);
snr_processed = snr(original, reconstructed - original);

y = snr_processed - snr_original;

end
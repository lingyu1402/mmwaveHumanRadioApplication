function [binIdx, zPeak, rd] = peak_track(raw, cfg, rlo, rhi)
% 逐 chirp 做距离 FFT，在 [rlo, rhi] 内找峰值，返回所在 bin 和复幅值
%
% binIdx 用 0 起算，距离 = binIdx * cfg.dR
% zPeak 是该 bin 上的复数样值，人走动时它的相位随距离连续变化，
% 静息状态下它的相位里就藏着呼吸和心跳

Ns = size(raw, 1);
M = size(raw, 2);
rd = fft(raw .* hamming(Ns), [], 1);

lo = max(1, floor(rlo / cfg.dR) + 1);
hi = min(Ns, ceil(rhi / cfg.dR) + 1);
if lo > hi
    error('距离搜索区间为空');
end

[~, ib] = max(abs(rd(lo:hi, :)), [], 1);
binIdx = (lo - 1) + ib - 1;
zPeak = rd(binIdx + 1 + (0:M-1)*Ns);

function raw = fmcw_echo(scat, cfg, noiseAmp)
% FMCW 中频回波生成
%
% scat(k).R    1 x M，第 k 个散射点在 M 个 chirp 时刻的瞬时距离(m)
% scat(k).amp  该散射点的幅度，缺省 1
% scat(k).phi  固定附加相位(rad)，缺省 0
% noiseAmp     复高斯噪声幅度，缺省 0
%
% 输出 raw: cfg.Ns x M 的快时间 x 慢时间矩阵
%
% 每个 chirp 内的距离变化在 100us 量级上远小于一个距离分辨率，
% 所以按"chirp 内距离冻结、chirp 间更新"处理。
% 剩余相位项 pi*S*tau^2 在 4m 内约 0.013rad，人体感知距离上可忽略；
% 到 10m 会涨到 0.2rad，做远距离高精度测相时要补回来。

if nargin < 3
    noiseAmp = 0;
end

M = numel(scat(1).R);
tfast = (0:cfg.Ns-1).' / cfg.fs;
raw = zeros(cfg.Ns, M);

for k = 1:numel(scat)
    Rk = scat(k).R(:).';
    if numel(Rk) ~= M
        error('scat(%d).R 长度与第一个散射点不一致', k);
    end
    amp = 1;
    phi0 = 0;
    if isfield(scat(k), 'amp'), amp = scat(k).amp; end
    if isfield(scat(k), 'phi'), phi0 = scat(k).phi; end

    fb = 2 * cfg.slope * Rk / cfg.c;
    phase = 2*pi * (tfast * fb) ...
          - 4*pi * (ones(cfg.Ns, 1) * Rk) / cfg.lam + phi0;
    raw = raw + amp * exp(1j * phase);
end

if noiseAmp > 0
    raw = raw + noiseAmp / sqrt(2) * ...
          (randn(cfg.Ns, M) + 1j * randn(cfg.Ns, M));
end

%% 项目6：微多普勒（Micro-Doppler）仿真
% 行走人体 / 小跑宠物 / 三叶风扇三种目标，比较 STFT 时频图上的微多普勒展宽
% 和慢时间包络的调制频率。人、宠物、风扇在这两个特征上分得开。
clear; close all; clc;

cfg = radar_cfg();
rng(23);
dur = 4.0;
M = round(dur / cfg.Tchirp);
t = (0:M-1) * cfg.Tchirp;
frf = cfg.frf;

sceneName = {'行走人体', '小跑宠物', '三叶风扇'};
res = repmat(struct(), 3, 1);

% STFT 窗长 0.3s：人体的摆动周期约 1s，0.3s 窗既能跟住时变又够分辨速度
winLen = round(0.3 * frf);
hop = round(0.04 * frf);
nfft = 4096;
fAxis = (-nfft/2:nfft/2-1) * frf / nfft;
vAxis = fAxis * cfg.lam / 2;

for s = 1:3
    scat = struct();
    switch s
        case 1
            % 步频 1Hz。腿摆幅 0.28m -> 峰值微多普勒 1.76m/s
            fg = 1.0;
            torso = 3.0 + 0.012*sin(2*pi*fg*t + 0.3);
            scat(1).R = torso;               scat(1).amp = 1.0;
            scat(2).R = 3.0 + 0.28*sin(2*pi*fg*t);        scat(2).amp = 0.55;
            scat(3).R = 3.0 - 0.28*sin(2*pi*fg*t);        scat(3).amp = 0.55;
            scat(4).R = 3.0 + 0.20*sin(2*pi*fg*t + 0.8);  scat(4).amp = 0.45;
            scat(5).R = 3.0 - 0.20*sin(2*pi*fg*t + 0.8);  scat(5).amp = 0.45;
            nScat = 5;
            rdTitle = '平均 R-D 谱：四肢摆动把单条多普勒线展成一片';

        case 2
            % 四条腿 3Hz 小跑，摆幅 0.09m，微多普勒峰值和人接近，
            % 但包络调制频率差三倍 —— 这是区分人和宠物的主要依据
            fg = 3.0;
            scat(1).R = 3.0 + 0.008*sin(2*pi*fg*t);       scat(1).amp = 0.55;
            ph = [0, pi, pi/2, -pi/2];
            for k = 1:4
                scat(k+1).R = 3.0 + 0.09*sin(2*pi*fg*t + ph(k));
                scat(k+1).amp = 0.35;
            end
            nScat = 5;
            rdTitle = '平均 R-D 谱：宠物四肢摆动把谱线展宽';

        case 3
            % 三叶旋转 4Hz。每片叶扫掠速度 4m/s，谱上是对称的两簇离散能量
            frot = 4.0;
            for k = 0:2
                scat(k+1).R = 3.0 + 0.16*cos(2*pi*frot*t + 2*pi*k/3);
                scat(k+1).amp = 0.5;
            end
            nScat = 3;
            rdTitle = '平均 R-D 谱：叶片旋转形成对称的两簇谱线';
    end

    scat(nScat+1).R = 4.4 + zeros(size(t));
    scat(nScat+1).amp = 0.22;

    raw = fmcw_echo(scat(1:nScat+1), cfg, 0.012);
    [~, ~, rd] = peak_track(raw, cfg, 2.3, 3.8);

    % 取目标所在的一段距离 bin。不同 bin 的 FFT 相位参考不共用，
    % 所以不能先相加再做 STFT；这里对每 bin 独立 STFT，然后在功率域累加
    binSel = round(2.4/cfg.dR) : round(3.7/cfg.dR);
    Z = rd(binSel, :);

    nWin = floor((M - winLen) / hop) + 1;
    S = zeros(nWin, nfft);
    w = hamming(winLen).';
    for i = 1:nWin
        idx = (i-1)*hop + (1:winLen);
        seg = fft(Z(:, idx) .* w, nfft, 2);
        S(i, :) = sum(abs(fftshift(seg, 2)).^2, 1);
    end
    tWin = ((0:nWin-1)*hop + winLen/2) * cfg.Tchirp;

    % 平均 R-D 图：长数据分成 256 chirp 一段，跨段做功率平均
    nSeg = floor(M / 256);
    rdAvg = zeros(cfg.Ns, 256);
    for i = 1:nSeg
        seg = rd(:, (i-1)*256 + (1:256));
        rdAvg = rdAvg + abs(fftshift(fft(seg .* hamming(256).', [], 2), 2)).^2;
    end
    rdAvg = rdAvg / nSeg;
    vRd = (-128:127) * cfg.lam / (2*256*cfg.Tchirp);

    % FFT 功率的绝对量级随点数走，画固定 caxis 会整片饱和。
    % 统一按各自峰值折成相对 dB，三个目标之间才可比较
    SdB = 10 * log10(S / max(S(:)) + 1e-12);
    rdB = 10 * log10(rdAvg / max(rdAvg(:)) + 1e-12);

    p = mean(S, 1);
    p0 = p;
    p0(abs(vAxis) < 0.15) = 0;              % 避开零速的墙面和支架回波
    vnz = vAxis(abs(vAxis) >= 0.15);
    pnz = p0(abs(vAxis) >= 0.15);
    vMean = sum(vnz .* pnz) / sum(pnz);
    spread = sqrt(sum((vnz-vMean).^2 .* pnz) / sum(pnz));
    [~, ip] = max(pnz);
    vPeak = abs(vnz(ip));

    % 慢时间包络的调制频率：四肢扫掠会让目标区总能量周期起伏
    env = sum(abs(Z), 1);
    env = env - mean(env);
    Ef = abs(fft(env .* hamming(M).')).^2;
    fEnv = (0:M/2-1) * frf / M;
    band = find(fEnv > 0.5 & fEnv < 15);
    [~, ifreq] = max(Ef(band));
    modHz = fEnv(band(ifreq));

    res(s).name = sceneName{s};
    res(s).rdTitle = rdTitle;
    res(s).SdB = SdB;
    res(s).tWin = tWin;
    res(s).vAxis = vAxis;
    res(s).p = p;
    res(s).rdB = rdB;
    res(s).vRd = vRd;
    res(s).spread = spread;
    res(s).vPeak = vPeak;
    res(s).modHz = modHz;
end

fprintf('\n目标       微多普勒展宽   峰值速度   包络调制频率\n');
fprintf('---------  ------------  --------  --------------\n');
for s = 1:3
    fprintf('%-9s  %11.3f  %7.3f  %12.2f\n', ...
        res(s).name, res(s).spread, res(s).vPeak, res(s).modHz);
end

for s = 1:3
    figure('Name', ['项目6 - ' res(s).name], 'Position', [100+25*s, 60+25*s, 1020, 850]);

    subplot(3, 1, 1);
    imagesc(res(s).tWin, res(s).vAxis, res(s).SdB.');
    axis xy; caxis([-45 0]); colormap(gca, parula); colorbar;
    ylim([-cfg.vmax cfg.vmax]);
    ylabel('速度 (m/s)');
    title([res(s).name ' - 微多普勒时频图（相对峰值 dB）']);

    subplot(3, 1, 2);
    imagesc(res(s).vRd, (0:cfg.Ns-1)*cfg.dR, res(s).rdB);
    axis xy; caxis([-40 0]); colormap(gca, parula); colorbar;
    xlim([-cfg.vmax cfg.vmax]); ylim([2.3 3.8]);
    ylabel('距离 (m)');
    title(res(s).rdTitle);

    subplot(3, 1, 3);
    plot(res(s).vAxis, 10*log10(res(s).p/max(res(s).p) + 1e-9), 'LineWidth', 1);
    grid on; xlim([-cfg.vmax cfg.vmax]); ylim([-50 0]);
    xlabel('速度 (m/s)'); ylabel('相对功率 (dB)');
    title(sprintf('速度谱：展宽 %.2f m/s，包络调制 %.2f Hz', ...
        res(s).spread, res(s).modHz));
end

%% 项目5：跌倒检测（径向速度路线 + 距离轨迹辅助）
% 坐下、跌倒、躺下、弯腰、快速弯腰五种轨迹。
% 主判据是速度峰值和速度衰减，距离轨迹只用来挡误报：
% 弯腰速度可以很快，但位移短而且会恢复，不能当跌倒。
clear; close all; clc;

cfg = radar_cfg();
rng(7);
Nchirp = 64;
nFrame = 250;
M = Nchirp * nFrame;
Tframe = Nchirp * cfg.Tchirp;
t = (0:M-1) * cfg.Tchirp;
step = @(x) (min(max(x, 0), 1)).^2 .* (3 - 2*min(max(x, 0), 1));

sceneName = {'坐下', '跌倒', '躺下', '弯腰', '快速弯腰'};
res = repmat(struct(), 5, 1);

for s = 1:5
    switch s
        case 1, R = 2.5 - 0.35*step(t/0.8);
        case 2, R = 2.5 - 1.40*step(t/0.45);
        case 3, R = 2.5 - 1.20*step(t/1.6);
        case 4, R = 2.5 - 0.50*step(t/0.5) + 0.50*step((t-0.6)/0.5);
        case 5, R = 2.5 - 0.50*step(t/0.33) + 0.50*step((t-0.45)/0.33);
    end

    scat(1).R = R;
    scat(1).amp = 1;
    scat(2).R = 4.4 + zeros(size(t));
    scat(2).amp = 0.22;

    raw = fmcw_echo(scat, cfg, 0.012);
    [binIdx, zPeak, ~] = peak_track(raw, cfg, 0.8, 3.2);
    Rf = double(binIdx(Nchirp:Nchirp:M)) * cfg.dR;

    % 帧间相位会绕回：跌倒峰值速度 4.7m/s 对应帧间相位 140rad，
    % 所以速度一律从 chirp 间差分取，vmax = lambda/(4*Tchirp) = 6.5m/s
    zExt = [zPeak, zPeak(end)];
    dphi = angle(zExt(2:end) .* conj(zExt(1:end-1)));
    vChirp = -dphi * cfg.lam / (4*pi*cfg.Tchirp);
    vF = mean(reshape(vChirp, Nchirp, nFrame), 1);

    tailN = round(0.5 / Tframe);
    vPeak = max(abs(vF));
    netDrop = Rf(1) - min(Rf);
    tailRms = rms(vF(end-tailN:end));
    recovery = Rf(end) - min(Rf);

    alarm = vPeak > 1.8 && netDrop > 0.6 && tailRms < 0.25 && recovery < 0.2;

    res(s).name = sceneName{s};
    res(s).Rf = Rf;
    res(s).vF = vF;
    res(s).tF = (0:nFrame-1) * Tframe;
    res(s).vPeak = vPeak;
    res(s).netDrop = netDrop;
    res(s).tailRms = tailRms;
    res(s).recovery = recovery;
    res(s).alarm = alarm;
    res(s).reason = '';
    if vPeak <= 1.8
        res(s).reason = '速度未达阈值';
    elseif netDrop <= 0.6
        res(s).reason = '位移不足';
    elseif tailRms >= 0.25 || recovery >= 0.2
        res(s).reason = '未进入静止 / 姿态已恢复';
    end
end

fprintf('\n场景       速度峰值   净位移   尾段速度RMS  姿态恢复  判决    抑制原因\n');
fprintf('---------  --------  -------  ------------  --------  ------  ----------------------\n');
for s = 1:5
    if res(s).alarm
        verdict = '报跌倒';
    else
        verdict = '不报';
    end
    fprintf('%-9s  %7.2f  %6.2f  %11.3f  %7.2f  %-6s  %s\n', ...
        res(s).name, res(s).vPeak, res(s).netDrop, res(s).tailRms, ...
        res(s).recovery, verdict, res(s).reason);
end

figure('Name', '项目5 - 跌倒检测时序', 'Position', [100, 80, 1150, 1080]);
for s = 1:5
    subplot(5, 2, 2*s-1);
    plot(res(s).tF, res(s).Rf, 'LineWidth', 1.2);
    set(gca, 'YDir', 'reverse');
    grid on; ylim([0.9 2.7]);
    ylabel('距离 (m)');
    title([res(s).name ' - 距离轨迹']);

    subplot(5, 2, 2*s);
    plot(res(s).tF, res(s).vF, 'LineWidth', 0.9); hold on;
    yline(1.8, 'r--'); yline(-1.8, 'r--');
    grid on; ylim([-6 6]);
    ylabel('径向速度 (m/s)');
    if res(s).alarm
        flag = '触发';
    else
        flag = '未触发';
    end
    title(sprintf('%s - %s', res(s).name, flag));
end
xlabel('时间 (s)');

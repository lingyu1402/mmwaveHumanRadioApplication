%% 项目4：人体微动 / 运动状态区分
% 静坐呼吸、翻身抬手、行走三种回波，比较帧间相位差、距离 bin 漂移和多普勒展宽。
% 判决结果直接决定后级用哪条处理链：微动走 unwrap 相位，大幅度动走差分/多普勒。
clear; close all; clc;

cfg = radar_cfg();
rng(19);
Nchirp = 64;
Tframe = Nchirp * cfg.Tchirp;
step = @(x) (min(max(x, 0), 1)).^2 .* (3 - 2*min(max(x, 0), 1));

sceneName = {'静坐呼吸', '翻身抬手', '行走'};
nFrameList = [700, 200, 200];
res = repmat(struct(), 3, 1);

for s = 1:3
    nFrame = nFrameList(s);
    M = Nchirp * nFrame;
    t = (0:M-1) * cfg.Tchirp;

    switch s
        case 1
            % 呼吸 3.5mm / 0.27Hz，心跳 0.6mm / 1.05Hz。
            % 两个幅度都远小于 7.8cm 的距离 bin，峰值 bin 不动，只能从相位里看
            Rbody = 2.5 + 3.5e-3*sin(2*pi*0.27*t) + 6e-4*sin(2*pi*1.05*t);
            scat(1).R = Rbody;
            scat(1).amp = 1;
            scat(2).R = 4.4 + zeros(size(t));   % 墙面静态回波
            scat(2).amp = 0.22;
            segT = 6.0;

        case 2
            % 0.55s 起抬手 0.42m，0.95s 起放回。峰值速度约 3.5m/s：
            % chirp 间相位增量约 2.1rad（差分还能用），帧间增量上百 rad（unwrap 必崩）
            Rbody = 2.5 + 2e-3*sin(2*pi*0.27*t);
            Rarm = 2.5 + 0.42*(step((t-0.55)/0.18) - step((t-0.95)/0.25));
            scat(1).R = Rbody;
            scat(1).amp = 1;
            scat(2).R = Rarm;
            scat(2).amp = 1.35;
            scat(3).R = 4.4 + zeros(size(t));
            scat(3).amp = 0.22;
            segT = 0.62;

        case 3
            % 朝雷达走 0.6m/s，1.9s 内距离从 2.5m 走到 1.35m，bin 单向漂移约 15 格
            Rbody = 2.5 - 0.6*t;
            Rleg1 = Rbody + 0.12*sin(2*pi*1.1*t);
            Rleg2 = Rbody - 0.12*sin(2*pi*1.1*t);
            scat(1).R = Rbody;
            scat(1).amp = 1;
            scat(2).R = Rleg1;
            scat(2).amp = 0.7;
            scat(3).R = Rleg2;
            scat(3).amp = 0.7;
            scat(4).R = 4.4 + zeros(size(t));
            scat(4).amp = 0.22;
            segT = 0.9;
    end

    raw = fmcw_echo(scat, cfg, 0.015);
    [binIdx, zPeak, rd] = peak_track(raw, cfg, 1.2, 4.0);

    % 帧末取样，得到和帧率一致的慢时间相位序列
    fi = Nchirp:Nchirp:M-1;
    zF = zPeak(fi);
    dphiF = angle(zF(2:end) .* conj(zF(1:end-1)));
    phiUnwrap = unwrap(angle(zF));

    % R-D 图取事件所在的一帧。FFT 前不做 MTI：静息微动一阶差分会被抹掉
    n0 = min(M-Nchirp+1, round(segT/cfg.Tchirp));
    rdSeg = rd(:, n0:n0+Nchirp-1);
    rdMap = fftshift(fft(rdSeg .* hamming(Nchirp).', [], 2), 2);
    vAxis = (-Nchirp/2:Nchirp/2-1) * cfg.lam / (2*Nchirp*cfg.Tchirp);

    bt = round(median(binIdx(n0:n0+Nchirp-1)));
    p = sum(abs(rdMap(max(1,bt-1):min(cfg.Ns,bt+1), :)).^2, 1);
    p = p(:).' - min(p);
    p(p < 0.08*max(p)) = 0;
    if sum(p) > 0
        vc = sum(vAxis.*p)/sum(p);
        dopSpread = sqrt(sum((vAxis-vc).^2.*p)/sum(p));
    else
        dopSpread = 0;
    end

    dbin = diff(double(binIdx));
    posCnt = sum(dbin > 0);
    negCnt = sum(dbin < 0);
    totalMove = posCnt + negCnt;
    backRatio = 0;
    if totalMove > 0
        backRatio = min(posCnt, negCnt) / totalMove;
    end

    res(s).name = sceneName{s};
    res(s).dur = M * cfg.Tchirp;
    res(s).phaseMax = max(abs(dphiF));
    res(s).phaseRms = rms(dphiF);
    res(s).phiStd = std(phiUnwrap - mean(phiUnwrap));
    res(s).binDrift = double(binIdx(end) - binIdx(1));
    res(s).binRange = double(max(binIdx) - min(binIdx));
    res(s).backRatio = backRatio;
    res(s).dopSpread = dopSpread;
    res(s).Rtrack = double(binIdx) * cfg.dR;
    res(s).phiF = phiUnwrap;
    res(s).dphiF = dphiF;
    res(s).tF = (0:numel(zF)-1) * Tframe;
    res(s).rdMap = rdMap;
    res(s).vAxis = vAxis;
    res(s).rPlot = (0:cfg.Ns-1)' * cfg.dR;

    % 判决：相位差分 RMS 区分微动和大幅度动。静息实测 0.15，
    % 运动在 0.7 以上，门限取 0.4 留了足够裕量。
    % 再用 bin 净漂移占总漂移范围的比例区分单向走动和短时往返
    driftRatio = abs(res(s).binDrift) / max(1, res(s).binRange);
    res(s).driftRatio = driftRatio;
    if res(s).phaseRms > 0.4
        if res(s).binDrift <= -4 && driftRatio > 0.7
            res(s).state = '行走';
            res(s).chain = '差分相位 / 多普勒跟踪';
        else
            res(s).state = '翻身或抬手';
            res(s).chain = 'chirp 间差分相位';
        end
    else
        res(s).state = '静坐微动';
        res(s).chain = 'unwrap 相位 + 谱分析';
    end
end

fprintf('\n场景        时长   相位峰值  相位RMS  bin漂移  单向度  多普勒展宽  判定          后级处理链\n');
fprintf('----------  -----  --------  -------  --------  ------  ----------  ------------  ----------------------\n');
for s = 1:3
    fprintf('%-10s  %.2f  %8.2f  %7.4f  %+7d  %6.2f  %9.3f  %-12s  %s\n', ...
        res(s).name, res(s).dur, res(s).phaseMax, res(s).phaseRms, ...
        res(s).binDrift, res(s).driftRatio, res(s).dopSpread, ...
        res(s).state, res(s).chain);
end

for s = 1:3
    figure('Name', ['项目4 - ' res(s).name], 'Position', [120+20*s, 100+20*s, 900, 640]);

    subplot(2, 1, 1);
    if strcmp(res(s).state, '静坐微动')
        plot(res(s).tF, res(s).phiF - mean(res(s).phiF), 'LineWidth', 1);
        ylabel('unwrap 相位 (rad)');
        title('慢时间相位：呼吸与心跳起伏');
    else
        plot(res(s).tF(1:end-1), res(s).dphiF, 'LineWidth', 0.8);
        yline(pi, 'r--'); yline(-pi, 'r--');
        ylabel('帧间相位差 (rad)');
        title('慢时间相位差分：超过 ±pi 后 unwrap 无法恢复');
    end
    grid on;
    xlabel('时间 (s)');
    legend(res(s).state, 'Location', 'best');

    subplot(2, 1, 2);
    imagesc(res(s).vAxis, res(s).rPlot, 20*log10(abs(res(s).rdMap)+1e-6));
    axis xy; caxis([-55 0]); colorbar;
    xlim([-cfg.vmax cfg.vmax]); ylim([1.2 4.0]);
    xlabel('径向速度 (m/s)'); ylabel('距离 (m)');
    title(sprintf('R-D 谱，目标附近展宽 %.3f m/s', res(s).dopSpread));
end

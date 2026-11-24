function cfg = radar_cfg()
% 77GHz FMCW 人体感知雷达参数
%
% 取 10 m 距离窗、理论距离分辨率 10 cm，反推 ADC 采样率和 chirp 斜率。
% 这套数值不是随意挑的：1.5GHz 带宽对应 c/(2B)=10cm，Tc=100us 时
% 10m 处的拍频约 1GHz，ADC 1MHz 采样在距离维上刚好覆盖到 10m。
% chirp 间隔留 50us 给 PLL 切换和 ADC 清空，帧长随使用的 chirp 数变化。

c = 3e8;
cfg.c = c;
cfg.f0 = 77e9;
cfg.lam = c / cfg.f0;

cfg.B = 1.5e9;
cfg.Tc = 100e-6;
cfg.Tidle = 50e-6;
cfg.Tchirp = cfg.Tc + cfg.Tidle;
cfg.slope = cfg.B / cfg.Tc;

cfg.fs = 1e6;
cfg.Ns = 128;

% 距离 bin 间隔比理论分辨率略小（7.8cm < 10cm），过采样一档，
% 这样峰值 bin 不会因为目标落在两个分辨率单元之间而跳来跳去
cfg.dR = c * cfg.fs / (2 * cfg.slope) / cfg.Ns;
cfg.Rmax = cfg.Ns * cfg.dR;

cfg.frf = 1 / cfg.Tchirp;
cfg.vmax = cfg.lam * cfg.frf / 4;

%% 标定加工好的NiTi管弹性模量E
% 拍摄NiTi管初始曲率图片，NiTi管插入一根SS棒后再次拍摄当前曲率图片.
% by Zhang Chao
% Date：2022/11/14

clear;clc;

%% NiTi管和SS管物理参数
Es = 193e9;         % SS材料弹性模量，Pa
OD_s = 1.93e-3;     % 所插入的SS棒外径

OD_n = 2.7e-3;      % NiTi管外径
ID_n = 2.4e-3;      % NiTi管内径
s = 175e-3;         % NiTi管实际弧长

%% 曲率计算 - 通过图片估算
% 初始曲率
alpha_1_fig = deg2rad(70.54);       % solidworks中测量所得
r1_fig = 1904.92;                   % solidworks中测量所得
s1_fig = r1_fig * alpha_1_fig;      % 图片计算弧长
ratio_1 = s/s1_fig;                 % 实际和图片中单位长度比例
r1_meas = ratio_1*r1_fig;
u_star_meas = 1/r1_meas

% 插入SS棒后的曲率(实际上，SS棒最大应变已达0.19%，超出弹性极限0.11%，所以只能用NiTi管插入来测试)
% 一个方法是测量抽出后的SS棒曲率变化，其应有塑性变形，可减去该部分效应
alpha_2_fig = deg2rad(19.66);       % solidworks中测量所得
r2_fig = 6825.34;                   % solidworks中测量所得
s2_fig = r2_fig * alpha_2_fig;      % 图片计算弧长
ratio_2 = s/s2_fig;                 % 实际和图片中单位长度比例
r2_meas = ratio_2*r2_fig;
u_meas = 1/r2_meas

%% 校准的弹性模量
En_calibrate = Es*OD_s^4/(OD_n^4-ID_n^4)*u_meas/(u_star_meas-u_meas)/1e9


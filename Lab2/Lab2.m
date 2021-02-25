clear;
clc;
close all;

% 1.1

load TRC_display.mat;
inputVals = 0:0.01:1;

colors = [0.80 0.00 0.00;
          0.00 0.80 0.00;
          0.00 0.00 0.80];
set(groot,'defaultAxesColorOrder',colors);

hold on;
plot(inputVals, TRCr);
plot(inputVals, TRCg);
plot(inputVals, TRCb);

% Too big of a difference between the channels will cause neutral colors to
% not actually look color neutral. For example with these TRCs neutral
% colors might look slightly blue/purple.

%% 1.2
clc;

load Ramp_display.mat;
load Ramp_linear.mat;

Ramp_display_compensated = compensateTRC(TRCr, TRCg, TRCb, Ramp_display);

imshow(Ramp_display_compensated);

%% 1.3
clc;

maxValue = max(max(max(Ramp_display)));

Ramp_display_gamma(:, :, 1) = maxValue * (Ramp_display(:, :, 1) / maxValue) .^ (1/2.1);
Ramp_display_gamma(:, :, 2) = maxValue * (Ramp_display(:, :, 2) / maxValue) .^ (1/2.4);
Ramp_display_gamma(:, :, 3) = maxValue * (Ramp_display(:, :, 3) / maxValue) .^ (1/1.8);

imshow(Ramp_display_gamma);

%% 2.1

load DLP.mat;

wavelengths = 400:5:700;

plot(wavelengths, DLP);

%% 2.2

load RGB_raw.mat;
load chips20.mat;
load illum.mat;
load xyz.mat;

Srgb_raw = DLP * RGB_raw; 

% Converting Srgb to CIEXYZ values.
k = 100 / (CIED65 * xyz(:, 2)); % normalization factor
XYZ_D65_ref = k * xyz' * (chips20 .* CIED65)';
Srgb_raw_XYZ_D65 = k * xyz' * Srgb_raw;

[mean1, max1] = deltaE(XYZ_D65_ref, Srgb_raw_XYZ_D65)

% The values are very large since we are using uncalibrated raw data.

%% 2.3

load RGB_cal.mat;

Srgb_cal = DLP * RGB_cal; 
Srgb_cal_XYZ_D65 = k * xyz' * Srgb_cal;

[mean2, max2] = deltaE(XYZ_D65_ref, Srgb_cal_XYZ_D65)

% With calibration the values are a bit better.

%% 3.1
% xyz is the observer and DLP the light from the projector.

ACRT = k * xyz' * DLP;

% ACRT is the relation between XYZ and RGB.

%% 3.2 

load XYZ_est.mat;

RGB_est = inv(ACRT) * XYZ_est;

SRGB = DLP * RGB_est;
XYZ_DLP_D65 = k * xyz' * SRGB;

[mean3, max3] = deltaE(XYZ_D65_ref, XYZ_DLP_D65)

% The values are very low when using the ACRT relation matrix.

%% 3.3

imshow(RGB_est);

% The RGB values include negative numbers but we dont have negative light.
% The method works better in theory than in real applications.

%% 3.4

norm_RGB_est = max(RGB_est, 0);
SRGB_norm = DLP * norm_RGB_est;
XYZ_DLP_D65_norm = k * xyz' * SRGB_norm;

[mean4, max4] = deltaE(XYZ_D65_ref, XYZ_DLP_D65_norm)

% The values are bigger now but more realistic.

%% 3.5 

plot_chrom_sRGB(ACRT);

% The gamut of the projector is smaller than that of sRGB standard and
% therefore colors that can't be produced have to be approximated to the
% closest color inside the gamut.

%% 3.6

plot(CIED65 .* chips20(1, :));
hold on;
plot(SRGB(:, 1));

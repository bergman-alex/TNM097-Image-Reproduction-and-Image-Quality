%% 1

% 1.1

clear;
clc;
close all;

load Ad.mat;
load Ad2.mat;
wavelengths = 400:5:700;

colors = [0.80 0.00 0.00;
          0.00 0.80 0.00;
          0.00 0.00 0.80];
set(groot,'defaultAxesColorOrder',colors)

plot(wavelengths, Ad);
figure;
plot(wavelengths, Ad2);

% The output from the two cameras wont necessarily be the same since the
% sensitivity functions differ, especially for higher wavelengths.

%%

% 1.2

clc;
close all;

load chips20.mat;
load illum.mat;

d1 = Ad' * (chips20(:, :) .* CIED65)'; % RGB output under outdoor light
showRGB(d1');

%save('RGB_raw_D65.mat', 'd1');

d2 = Ad2' * (chips20(:, :) .* CIED65)'; % RGB output under outdoor light
showRGB(d2');

% We can see a clear difference. For example, Ad2 is more
% sensitive to green than Ad, which is apparent.

%% 2

% 2.1

clc;
close all;

load Ad.mat;
load Ad2.mat;
load chips20.mat;
load illum.mat;
load RGB_raw_D65.mat;
load RGB_cal_D65.mat;

d1_n = ones(3, 1) ./ (Ad'  * ones(1, 61)'); % normalization factors
d2_n = ones(3, 1) ./ (Ad2' * ones(1, 61)');

% The red values of the normalization factors are higher since the red
% sensitivity is the lowest.

%%

% 2.2

clc;
close all;

dPrime1 = d1 .* d1_n; % calibrated camera responses
dPrime2 = d2 .* d2_n;

%save('RGB_cal_D65.mat', 'dPrime1');

showRGB(dPrime1');
showRGB(dPrime2');

% The two signals are now very similar

%%

% 2.3

clc;
close all;

wavelengths = 400:5:700;

plot(wavelengths, CIEA);
figure;
plot(wavelengths, CIED65);

% One difference seems to be that the indoor light has more warm colors
% such as yellow and red, while the outdoor light has bluer colors.
% The indoor light also has a much smoother curve which i would guess is
% because it's man-made  while the outdoor light is sampled from the sun.

%%

% 2.4

clc;
close all;

load RGB_cal_D65.mat;

d1_indoor = Ad' * (chips20(:, :) .* CIEA)'; % RGB output under indoor light
dPrime_indoor = d1_indoor .* d1_n; % calibrated camera response

showRGB(dPrime1'); % calibrated outdoor
showRGB(dPrime_indoor'); % calibrated indoor

% The difference in illumination source has a big impact on the result.

%%

% 2.5

clc;
close all;

d1_indoor_n  = ones(3, 1) ./ (Ad' * (ones(1, 61) .* CIEA)'); % normalization factors with light calibration included
d1_outdoor_n = ones(3, 1) ./ (Ad' * (ones(1, 61) .* CIED65)');

d1_indoor = Ad' * (chips20(:, :) .* CIEA)'; % RGB output under indoor light
d1_outdoor = Ad' * (chips20(:, :) .* CIED65)'; % RGB output under indoor light

dPrime_indoor_cal = d1_indoor .* d1_indoor_n; % response calibrated according to light
dPrime_outdoor_cal = d1_outdoor .* d1_outdoor_n;

showRGB(dPrime_indoor_cal');
showRGB(dPrime_outdoor_cal');

%% 3.1

clear;
clc;
close all;

load chips20.mat;
load illum.mat;
load xyz.mat;
load Ad.mat;
load M_XYZ2RGB.mat;
load RGB_cal_D65.mat;
dPrime_rgb_d65 = dPrime1';
clear dPrime1;

d = xyz' * (chips20(:, :) .* CIED65)';
d_n = 100 ./ (CIED65*xyz(:,2)); % normalization factor

dPrime_xyz_d65 = (d .* d_n)';

save('XYZ_D65_ref.mat', 'dPrime_xyz_d65');

%% 3.2

xyz_est = (inv(M_XYZ2RGB) * dPrime_rgb_d65')'; % converting to XYZ using the inverse matrix

[est_L, est_a, est_b] = xyz2lab(xyz_est(:,1), xyz_est(:,2), xyz_est(:,3)); % get CIELab values of the estimate
[ref_L, ref_a, ref_b] = xyz2lab(dPrime_xyz_d65(:,1), dPrime_xyz_d65(:,2), dPrime_xyz_d65(:,3)); % get CIELab values of the reference

deltaE = sqrt((est_L - ref_L).^2 + (est_a - ref_a).^2 + (est_b - ref_b).^2); % calculate euclidean distance

maxDiff = max(deltaE);
meanDiff = mean(deltaE);

%% 3.3

plot(400:5:700, Ad);
figure;
plot(400:5:700, xyz);

%% 3.4

A = pinv(dPrime_rgb_d65) * dPrime_xyz_d65;

xyz_est2 = (dPrime_rgb_d65 * A);

[est2_L, est2_a, est2_b] = xyz2lab(xyz_est2(:,1), xyz_est2(:,2), xyz_est2(:,3)); % get CIELab values of the estimate

deltaE2 = sqrt((est2_L - ref_L).^2 + (est2_a - ref_a).^2 + (est2_b - ref_b).^2); % calculate euclidean distance

maxDiff = max(deltaE2);
meanDiff = mean(deltaE2);

% The difference is now much smaller with the empirical technique, since
% the earlier estimation was not colorimetric.

%% 3.5 

A_optimal = Optimize_poly(dPrime_rgb_d65', dPrime_xyz_d65');
xyz_est3 = Polynomial_regression(dPrime_rgb_d65', A_optimal)';

[est3_L, est3_a, est3_b] = xyz2lab(xyz_est3(:,1), xyz_est3(:,2), xyz_est3(:,3)); % get CIELab values of the estimate

deltaE3 = sqrt((est3_L - ref_L).^2 + (est3_a - ref_a).^2 + (est3_b - ref_b).^2); % calculate euclidean distance

maxDiff = max(deltaE3);
meanDiff = mean(deltaE3);

% The difference is even smaller when using polynomial regression since
% it's optimized for this specific scenario.







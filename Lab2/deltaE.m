function [mean_E,max_E]=deltaE(XYZ_ref,XYZ_est)

%calculate the CIELAB
[L1,a1,b1] = xyz2lab(XYZ_ref(1,:),XYZ_ref(2,:),XYZ_ref(3,:));
[L2,a2,b2] = xyz2lab(XYZ_est(1,:),XYZ_est(2,:),XYZ_est(3,:));

deltaE = sqrt((L1-L2).^2 + (a1-a2).^2 + (b1-b2).^2);
mean_E = mean(deltaE);
max_E = max(deltaE);

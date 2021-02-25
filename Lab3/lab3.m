%% Color Gamut
load('Dell.mat');
load('Inkjet.mat')

plot_chrom(XYZdell, 'b');
hold on 
plot_chrom(XYZinkjet, 'g');
xlabel('red');
ylabel('green');

% % The red line is the whole colorspace
% The dell should give better result in showing red

%% 2. Metrics
% 2. Mathemacial metrics
% 2.1 Grayscale image (MSE/SNR)
% 2.1.1 Interpolation

% Read image
image = im2double(imread('peppers_gray.tif'));

% Nearest
nearest= imresize(imresize(image,0.25,'nearest'),4,'nearest');
imshow(nearest);
figure

% Bilinear
bilinear= imresize(imresize(image,0.25,'bilinear'),4,'bilinear');
imshow(bilinear);
figure

% Bicubic
bicubic= imresize(imresize(image,0.25,'bicubic'),4,'bicubic');
imshow(bicubic);

% SNR: Signal to Noise Ratio
SNR_neareset=snr(image,image-nearest)
SNR_bilinear=snr(image,image-bilinear)
SNR_bicubic=snr(image,image-bicubic)

% MSE: Mean Square Error 
MSE_neareset=immse(image,nearest)
MSE_bilinear=immse(image,bilinear)
MSE_bicubic=immse(image,bicubic)

% These metrics are very fast and simple to implement but
% usually not well correlated with perceived image quality.

% The bicubic is the best reproduction when it comes to my judgement

%% 2.1.2 Halftoning
image = im2double(imread('peppers_gray.tif'));

% Halftone the image with error diffusion
dither_img = dither(image);
level = graythresh(image);
thresh_img = image >= level;

% Convert to double
dither_img = double(dither_img);
thresh_img = double(thresh_img);

% SNR
SNR_thresh=snr(image,image - thresh_img)
SNR_dither=snr(image,image - dither_img)

% MSE
MSE_thresh = immse(image,thresh_img)
MSE_dither = immse(image,dither_img)

images = cat(2,thresh_img,dither_img);
montage(images);


%% 2.2 Color Image Delta Eab
imageRGB = im2double(imread('peppers_color.tif'));

% Create referens Lab image
imLab = rgb2lab(imageRGB);

% Create threshold
% Halftone the color image by thresholding it with 0.5
ThreshRGB(:,:,1) = imageRGB(:,:,1)>=0.5;
ThreshRGB(:,:,2) = imageRGB(:,:,2)>=0.5;
ThreshRGB(:,:,3) = imageRGB(:,:,3)>=0.5;

% Create dither image
DitherRGB(:,:,1) = dither(imageRGB(:,:,1));
DitherRGB(:,:,2) = dither(imageRGB(:,:,2));
DitherRGB(:,:,3) = dither(imageRGB(:,:,3));

% Convert to double
DitherRGB = double(DitherRGB);
ThreshRGB = double(ThreshRGB);

% Convert to Lab
DitherLAB =rgb2lab(DitherRGB);
ThreshLAB =rgb2lab(ThreshRGB);

% Calculate the difference
% A pixel-by-pixel difference in CIELab
deltaETresh= 0;
for i = 1:1:512
    for j = 1:1:512
    
    deltaETresh = deltaETresh + sqrt((ThreshLAB(i,j,1)-imLab(i,j,1))^2 +  (ThreshLAB(i,j,2)-imLab(i,j,2))^2+ (ThreshLAB(i,j,3)-imLab(i,j,3))^2);
    
    end

end
deltaETresh = deltaETresh/(512*512);


deltaEDith= 0;
for i = 1:1:512
    for j = 1:1:512
    
    deltaEDith = deltaEDith + sqrt((DitherLAB(i,j,1)-imLab(i,j,1))^2 +  (DitherLAB(i,j,2)-imLab(i,j,2))^2 + (DitherLAB(i,j,3)-imLab(i,j,3))^2);
    
    end

end
deltaEDith = deltaEDith/(512*512);

% Show images 
figure(1)
imshow(ThreshRGB)
figure(2)
imshow(DitherRGB)

deltaETresh
deltaEDith

%% 3  Mathematical metrics involving HVS
img = im2double(imread('peppers_gray.tif'));

imgThres = im2double(img >= 0.5);
imgDith = im2double(dither(img));

snrDith = snr_filter(img, img-imgDith);
snrTresh = snr_filter(img, img-imgThres);

imshow(img)
figure
imshow(imgThres)
figure
imshow(imgDith)

% Metrics 
snrDith
snrTresh

%% 3.1

imageRGB = im2double(imread('peppers_color.tif'));

% Convert to double
DitherRGB = double(DitherRGB);
ThreshRGB = double(ThreshRGB);

% Filter of what an eye sees
eye_filter=MFTsp(15,0.0847,500);

% Referens
ref(:,:,1) =conv2(imageRGB(:,:,1),eye_filter,'same');
ref(:,:,2) =conv2(imageRGB(:,:,2),eye_filter,'same');
ref(:,:,3) =conv2(imageRGB(:,:,3),eye_filter,'same');
ref(:,:,1) = (ref(:,:,1)>0).*ref(:,:,1);
ref(:,:,2) = (ref(:,:,2)>0).*ref(:,:,2);
ref(:,:,3) = (ref(:,:,3)>0).*ref(:,:,3);
ref = rgb2lab(ref);

% Dither
DitherHVS(:,:,1) =conv2(DitherRGB(:,:,1),eye_filter,'same');
DitherHVS(:,:,2) =conv2(DitherRGB(:,:,2),eye_filter,'same');
DitherHVS(:,:,3) =conv2(DitherRGB(:,:,3),eye_filter,'same');

DitherHVS(:,:,1) = (DitherHVS(:,:,1)>0).*DitherHVS(:,:,1);
DitherHVS(:,:,2) = (DitherHVS(:,:,2)>0).*DitherHVS(:,:,2);
DitherHVS(:,:,3) = (DitherHVS(:,:,3)>0).*DitherHVS(:,:,3);
DitherHVS = rgb2lab(DitherHVS);

% Halftone
HalftoneHVS(:,:,1) =conv2(ThreshRGB(:,:,1),eye_filter,'same');
HalftoneHVS(:,:,2) =conv2(ThreshRGB(:,:,2),eye_filter,'same');
HalftoneHVS(:,:,3) =conv2(ThreshRGB(:,:,3),eye_filter,'same');

HalftoneHVS(:,:,1) = (HalftoneHVS(:,:,1)>0).*HalftoneHVS(:,:,1);
HalftoneHVS(:,:,2) = (HalftoneHVS(:,:,2)>0).*HalftoneHVS(:,:,2);
HalftoneHVS(:,:,3) = (HalftoneHVS(:,:,3)>0).*HalftoneHVS(:,:,3);
HalftoneHVS = rgb2lab(HalftoneHVS);

% Delta e for dither
DeltaEdither = sqrt((DitherHVS(:,:,1)-ref(:,:,1)).^2+(DitherHVS(:,:,2)-ref(:,:,2)).^2+(DitherHVS(:,:,3)-ref(:,:,3)).^2);
DeltaEditherMetric = sum(sum(DeltaEdither))/(size(DeltaEdither,1)*size(DeltaEdither,2));

% Delta e for thres
DeltaEThresh = sqrt((HalftoneHVS(:,:,1)-ref(:,:,1)).^2+(HalftoneHVS(:,:,2)-ref(:,:,2)).^2+(HalftoneHVS(:,:,3)-ref(:,:,3)).^2);
DeltaEThreshMetric = sum(sum(DeltaEThresh))/(size(DeltaEThresh,1)*size(DeltaEThresh,2));

imshow(ref)
figure
imshow(HalftoneHVS)
figure
imshow(DitherHVS)

% Metrics
DeltaEditherMetric
DeltaEThreshMetric

%% 4.1 S-CIELab as a full-reference metric

%Load image
img = imread('peppers_color.tif');

imgNe = imresize(imresize(img,0.25,'nearest'),4,'nearest');
imgLi = imresize(imresize(img,0.25,'bilinear'),4,'bilinear');
imgBi = imresize(imresize(img,0.25,'bicubic'),4,'bicubic');

%Convert to XYZ
imgXYZ = rgb2xyz(img);
nearestXYZ = rgb2xyz(imgNe);
linearXYZ = rgb2xyz(imgLi);
bicubicXYZ = rgb2xyz(imgBi);

%Scielab
%Sample
% 100/2.5 = 1 m
sampDegree = visualAngle(-1,100/2.5,72,1);
whitePoint = [95.05,100,108.9];

near = scielab(sampDegree,imgXYZ,nearestXYZ,whitePoint,'xyz');
linear = scielab(sampDegree,imgXYZ,linearXYZ,whitePoint,'xyz');
bicubic = scielab(sampDegree,imgXYZ,bicubicXYZ,whitePoint,'xyz');

nearMean = mean(mean(near));
linearMean = mean(mean(linear));
bicubicMean = mean(mean(bicubic));

imshow(img);
figure
imshow(imgNe);
figure
imshow(imgLi);
figure
imshow(imgBi);

% Metrics
nearMean
linearMean
bicubicMean

%% 4.2 S-CIELab as a no-reference metric
% We want to test this metric to evaluate the graininess of a 
% number of color halftones.

% In c1, the cyan and magenta channels have been halftoned independently
% in c2 the channels have been halftoned dependently.

% colorhalftones have five color halftone patches, called c1 through c5.
load('colorhalftones.mat');
sampDegree = visualAngle(-1,500/2.5,300,1);

c1xyz = rgb2xyz(c1);
c2xyz = rgb2xyz(c2);

c1Lab = scielab(sampDegree,c1xyz);
c2Lab = scielab(sampDegree,c2xyz);

c1Sum = sum(std(c1Lab(:,:,1))) + sum(std(c1Lab(:,:,2)))+sum(std(c1Lab(:,:,3)));
c2Sum = sum(std(c2Lab(:,:,1)))+sum(std(c2Lab(:,:,2)))+sum(std(c2Lab(:,:,3)));

imshow(c1);
figure;
imshow(c2);

% Metrics
c1Sum
c2Sum

%% 4.2.2 

load('colorhalftones.mat');
sampDegree = visualAngle(-1,500/2.5,300,1);

c3xyz = rgb2xyz(c3);
c4xyz = rgb2xyz(c4);
c5xyz = rgb2xyz(c5);

c3Lab = scielab(sampDegree,c3xyz);
c4Lab = scielab(sampDegree,c4xyz);
c5Lab = scielab(sampDegree,c5xyz);

c3Sum = sum(std(c3Lab(:,:,1))) + sum(std(c3Lab(:,:,2)))+sum(std(c3Lab(:,:,3)));
c4Sum = sum(std(c4Lab(:,:,1)))+sum(std(c4Lab(:,:,2)))+sum(std(c4Lab(:,:,3)));
c5Sum = sum(std(c5Lab(:,:,1)))+sum(std(c5Lab(:,:,2)))+sum(std(c5Lab(:,:,3)));

imshow(c3);
figure;
imshow(c4);
figure;
imshow(c5);

% Metrics
c3Sum
c4Sum
c5Sum

%% 5 SSIM structural similarity
% Part A

img = im2double(imread('peppers_gray.tif'));
imgA1 = img;
imgA2 = img;

% Distortion 1: add +0.1 to odd rows and -0.1 to the other rows
% in the original image.
for i= 1:1:512
   for j = 1:1:512
    if( mod(i,2) == 0)
        imgA1(i,j) = img(i,j) - 0.1;
    else
        imgA1(i,j) = img(i,j) + 0.1;
    end
   end
end

% Distortion 2: add 0.1 to the upper half of the original image and 
% subtract 0.1 from the lower half of it.
for i= 1:1:512
   for j = 1:1:512
    if( i< 256)
        imgA2(i,j) = img(i,j) + 0.1;
    else
        imgA2(i,j) = img(i,j) - 0.1;
    end
   end
end

% Both distortions will exactly give the same
% SNR when compared to the original image (check it), why?
SNR_A1=snr(img,img-imgA1);
SNR_A2=snr(img,img-imgA2);

% Use ssim to check how similar these two distorted images are to the original.
[A1, SSIMA1] = ssim(img, imgA1);
[A2, SSIMA2] = ssim(img, imgA2);

imshow(SSIMA1);
figure
imshow(SSIMA2);

% Metrics
A1
A2

%% 5 SSIM structural similarity
% Part B
img = im2double(imread('peppers_gray.tif'));

% Distortion 1: Add noise to the original image
dist1=img+0.2*(rand(size(img))-0.5);

% Distortion 2: Use a Gaussian filter on the original image
f=fspecial('gauss',21,10);
dist2 = imfilter(img,f);

% Get snr values
SNR_B1=snr(img,img-dist1);
SNR_B2=snr(img,img-dist2);

[B1,SSIMB1] = ssim(img, dist1);
[B2,SSIMB2] = ssim(img, dist2);

% Show distortion maps
imshow(img);
figure
imshow(SSIMB1);
figure
imshow(SSIMB2);

% Metrics
B1
B2

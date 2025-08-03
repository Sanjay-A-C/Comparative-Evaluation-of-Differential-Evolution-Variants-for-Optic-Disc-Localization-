% Differential Evolution Algorithm for localization of Optic Disc in
% Retinal images
warning("off","all");
clear all;
close all;
clc;
tic;
ngen = 200; % Number of generations
npop = 20; % Number of chromosomes in population
im = imread('/MATLAB Drive/diaretdb1/diaretdb1_image038.png'); 

% Display the original image
%figure;
%imshow(im);
%title('Original Image');

x=rgb2gray(im);
[m n]=size(x);
x = medfilt2(x,[32 32]);
x = adapthisteq(x);
x1 = medfilt2(x, [110 110]); %previously it was 110 110
x=x-x1;
x = medfilt2(x, [32 32]); %previous it was 12 12
%average filtering
C1=fspecial('average',[40 40]);%previous it was 20 20
d101=imfilter(x,C1);
x=d101;

%figure;
%imshow(x);
%title('Preprocessed Image');

[pop(:,:,1),theta] = QDEsolset_sir(npop,m,n); % Generation of Population

unfitinit = zeros(npop, 1);
for j = 1:npop
[unfitinit(j,1)] = QDEfitnessvalue1entropy_sir(x, pop(j,:)); % Finding fitness value for each chromosome (initial population)
end

a=1:npop;
% SHADE Parameters
H=20;
M_F = 0.5 * ones(1, H);
M_Cr = 0.5 * ones(1, H);
alpha =0.3;


% Differential Evolution Algorithm starts
qtableunfitval1(1:npop, 1:ngen) = zeros(npop, ngen);
tic;
for i = 1:ngen
if(i==1)
k=1;
else
k=i-1;
end

% Generate F and Cr using historical memory
    r = randi(H, [npop, 1]);
    F = M_F(r) + 0.1 * randn(npop, 1);
    F(F > 1) = 1; F(F < 0) = 0;
    Cr = M_Cr(r) + 0.1 * randn(npop, 1);
    Cr(Cr > 1) = 1; Cr(Cr < 0) = 0;

for j=1:npop
b=randperm(npop);
c=find(b~=j);
r1=a(b(c(1)));
r2=a(b(c(2)));
r3=a(b(c(3)));
%dv = round(abs(pop(r1,:,k) + (F(j) .* (-(pop(r2,:,k) - pop(r3,:,k))))));
dv = QDE_mutation(pop(:,:,k), F(j), r1, r2, r3);

uv = QDEopticdecross_sir(pop(j,:,k),dv,Cr(j),m,n);

    if i > 10
        uvv = QDElocalsearch(uv, m, n, i, ngen);
    else
        uvv = uv;
    end

q1=QDEfitnessvalue1entropy_sir(x,uvv);
q2=QDEfitnessvalue1entropy_sir(x,pop(j,:,k));

if(q1>=q2)
pop(j,:,i)= round(uvv(:,:));
qtableunfitval1(j,i)=q1;
else
pop(j,:,i) = round(pop(j,:,k));
qtableunfitval1(j,i)=q2;
end

end

[qtableunfitval2(:,i),qtableindex(:,i)]=sort(qtableunfitval1(:,i),1,'descend');
% Update Historical Memory
    success = find(qtableunfitval1(:,i) > qtableunfitval1(:, max(i-1,1)));
    if ~isempty(success)
        k = mod(k-1, H) + 1;  % Ensure k cycles within H
        M_F(k) = mean(F(success));  % Replace only index k
        M_Cr(k) = mean(Cr(success));  % Replace only index k
    end
%disp(['Generation: ', num2str(i), ', H: ', num2str(H)]);
%disp(['Size of M_F: ', num2str(size(M_F))]);
%disp(['Size of M_Cr: ', num2str(size(M_Cr))]);

i;
if i == ngen
break;
end
end
% Differential Evolution Algorithm ends
t=toc;

% The following part is used to find best chromosome and its fitness value
[unfitsort(1:npop, i), index(:,i)]=sort(qtableunfitval1(:,i),1,'descend');
qtsort(1:npop,:,i)=pop(index(:,i),:,i);

qtunfitinitsort = [];
unfitinitsort = [];
npop = [];
npop = numel(unfitsort(:,i));
qtunfitinitsort(1:npop,:) = qtsort(:,:,i);
unfitinitsort(1:npop, i) = unfitsort(:,i);

unfitinitsort2 = unfitinitsort(:,i);
qtunfitinitsort2 = qtunfitinitsort;

averageunfitval1(1,:)=mean(qtableunfitval1);
bestunfitval(1,:)=qtableunfitval2(npop,:);
b1=bestunfitval(1,:);
% end of part
% This part is used find the change between consecutive generations
leap=0;
for i=2:length(b1)
if (b1(i-1)<b1(i))
leap=leap+1;
end
if (i==50)
leap50=leap;
elseif (i==100)
leap100=leap;
elseif(i==150)
    leap150=leap;
elseif(i==200)
    leap200=leap;
end
end
%end of part
optimqtable1 = qtunfitinitsort2(end,:); % best chromosome. This gives lefttop of optic disc. This is calculated through DE algorithm
%disp(['Leap Value: ', num2str(leap)]);

groundtruth=[1040 530 10 10]; % left top of optic disc. [xcordinates ycordinates r1 r2]. This is calculated manually for the input image
distance=sqrt(((groundtruth(1,1)-optimqtable1(1,1))^2)+((groundtruth(1,2)- optimqtable1(1,2))^2)); % Distance between ground truth and OD using DE
disp('distance : ');
disp(distance);
%results(1,:)=[optimqtable1(1,1);optimqtable1(1,2);optimqtable1(1,3);optimqtable1(1,4);leap25; leap50; t;distance];
[x1 y1 width height] = squaretable(optimqtable1);

figure();
imshow(im);
title('Optic Disc Localisation');
rectangle('Position',[x1,y1,width,height],'Curvature',[0.8,0.4],'LineWidth',2 , 'LineStyle','--');

x1 = round(abs(optimqtable1(1,1)));
y1 = round(abs(optimqtable1(1,2)));
width = round(abs(optimqtable1(1,3)));
height = round(abs(optimqtable1(1,4)));
% Validate width and height
if width > 0 && height > 0
    detected_region = im(y1:y1+height, x1:x1+width);
    
    % Compute metrics
    mean_intensity = mean(detected_region(:));
    std_dev = std(double(detected_region(:)));
else
    mean_intensity = NaN;
    std_dev = NaN;
    fprintf('Warning: Detected region is invalid!\n');
end

detected_area = width * height;
% Circularity
perimeter = 2 * (width + height); % Approximate perimeter
circularity = (4 * pi * detected_area) / (perimeter^2);

results(1,:)=[optimqtable1(1,1);optimqtable1(1,2);optimqtable1(1,3);optimqtable1(1,4);0; distance;0; mean_intensity;std_dev;circularity;leap50; leap100; leap150; leap200];

% Display Results
fprintf('Mean Intensity: %f\n', mean_intensity);
fprintf('Standard Deviation: %f\n', std_dev);
fprintf('Circularity: %f\n', circularity);
fprintf('%d\t , %d\t,%d\t,%d\t',x1,y1,width,height);
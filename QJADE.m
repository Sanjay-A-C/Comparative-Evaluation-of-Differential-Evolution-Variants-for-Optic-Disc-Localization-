% Differential Evolution Algorithm for localization of Optic Disc in
% Retinal images
warning("off","all");
clear all;
close all;
clc;
tic;
ngen = 200; % Number of generations
npop = 20; % Number of chromosomes in population
im = imread('/MATLAB Drive/diaretdb1/diaretdb1_image084.png'); 

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
% JADE Parameters
mu_F = 0.5;
mu_CR = 0.5;
archive = zeros(0, size(pop, 2));
archive_size = npop; % Maximum size of archive

% Differential Evolution Algorithm starts
qtableunfitval1(1:npop, 1:ngen) = zeros(npop, ngen);

tic;
for i = 1:ngen
if(i==1)
k=1;
else
k=i-1;
end

F = mu_F + 0.1 * randn(npop, 1);
CR = mu_CR + 0.1 * randn(npop, 1);
F = max(0, min(F, 1));
CR = max(0, min(CR, 1));

for j=1:npop
b=randperm(npop);
c=find(b~=j);
%r1=a(b(c(1)));
%r2=a(b(c(2)));
%r3=a(b(c(3)));
%dv=round(abs(pop(r1,:,k)+F*(-(pop(r2,:,k)-pop(r3,:,k)))));

theta = rand() * pi;
Q_theta = sin(theta).^2;

idx = randperm(npop);
        idx(idx == j) = [];
        r1 = idx(1); r2 = idx(2); r3 = idx(3);
        
        if ~isempty(archive)
            r4 = randi(size(archive,1));
            dv = round(abs(pop(r1,:,k) + F(j) * (pop(r2,:,k) - archive(r4,:)) + F(j) * (pop(r3,:,k) - pop(r1,:,k))+Q_theta));
        else
            dv = round(abs(pop(r1,:,k) + F(j) * (pop(r2,:,k) - pop(r3,:,k))+Q_theta));
        end
        
% Check and correct out-of-bound values:
    if dv(1) < 150 || dv(1) > 1350
        dv(1) = pop(r1, 1);  % Take from current population if out of bounds
    end
    if dv(2) < 150 || dv(2) > 1002
        dv(2) = pop(r1, 2);  % Take from current population if out of bounds
    end
    if dv(3) < -10 || dv(3) > 10
        dv(3) = pop(r1, 3);  % Take from current population if out of bounds
    end
    if dv(4) < -10 || dv(4) > 10
        dv(4) = pop(r1, 4);  % Take from current population if out of bounds
    end

        uv = QDEopticdecross_sir(pop(j,:,k), round(abs(dv)), CR(j), m, n);
%uv = DEopticdecross_sir(pop(j,:,k),dv,cro,m,n);

    if i > 10
        uvv = QDElocalsearch(uv, m, n, i, ngen);
    else
        uvv = uv;
    end

q1=QDEfitnessvalue1entropy_sir(x,uvv);
q2=QDEfitnessvalue1entropy_sir(x,pop(j,:,k));

% Selection and Archive Update
        if q1 >= q2
            % If new solution is better, store the old one in archive
            if size(archive, 1) < archive_size
                archive = [archive; pop(j,:,k)];
            else
                archive(randi(archive_size), :) = pop(j,:,k);
            end
            pop(j,:,i) = round(uvv(:,:));
            qtableunfitval1(j,i) = q1;
        else
            pop(j,:,i) = round(pop(j,:,k));
            qtableunfitval1(j,i) = q2;
        end

end

[qtableunfitval2(:,i),qtableindex(:,i)]=sort(qtableunfitval1(:,i),1,'descend');
% Update Adaptive Parameters
    if sum(qtableunfitval1(:,i) > qtableunfitval1(:,k)) > 0
        mu_F = (1 - 0.1) * mu_F + 0.1 * mean(F(qtableunfitval1(:,i) > qtableunfitval1(:,k)));
        mu_CR = (1 - 0.1) * mu_CR + 0.1 * mean(CR(qtableunfitval1(:,i) > qtableunfitval1(:,k)));
    end

    % Debugging Statements
    %disp(['Generation: ', num2str(i)]);
    %disp(['Archive Size: ', num2str(size(archive, 1))]);
    %disp(['mu_F: ', num2str(mu_F), ' mu_CR: ', num2str(mu_CR)]);
    %disp(['Generation: ', num2str(i)]);
    %disp(['Updated mu_F: ', num2str(mu_F), ' Updated mu_CR: ', num2str(mu_CR)]);
    %disp(['First 5 F values: ', num2str(F(1:5)')]);
    %disp(['First 5 CR values: ', num2str(CR(1:5)')]);

i;
if i == ngen
break;
end
end
% Differential Evolution Algorithm ends
t=toc;
%disp(archive);
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

groundtruth=[1180 560 10 10]; % left top of optic disc. [xcordinates ycordinates r1 r2]. This is calculated manually for the input image
distance=sqrt(((groundtruth(1,1)-optimqtable1(1,1))^2)+((groundtruth(1,2)- optimqtable1(1,2))^2)); % Distance between ground truth and OD using DE
disp('distance : ');
disp(distance);

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

results = [optimqtable1(1,1),optimqtable1(1,2),optimqtable1(1,3),optimqtable1(1,4),0,distance,0,mean_intensity, std_dev, circularity,leap50, leap100,leap150,leap200];
% Display Results
fprintf('Mean Intensity: %f\n', mean_intensity);
fprintf('Standard Deviation: %f\n', std_dev);
fprintf('Circularity: %f\n', circularity);
fprintf('%d\t , %d\t,%d\t,%d\t',x1,y1,width,height);

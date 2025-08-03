% Differential Evolution Algorithm with Quantum Concept for localization of Optic Disc in
% Retinal images
warning("off","all");
clear all;
close all;
clc;
tic;
 
ngen = 200;  % Number of generations
npop = 20; % Number of chromosomes in population
 
im=imread('/MATLAB Drive/diaretdb1/diaretdb1_image049.png'); % read the image
 
x=rgb2gray(im);
[m n]=size(x);
x = medfilt2(x,[32 32]);
x = adapthisteq(x);
 
x1 = medfilt2(x, [110 110]);
x=x-x1;
 
x = medfilt2(x, [32 32]);
%average filtering
C1=fspecial('average',[40 40]);
d101=imfilter(x,C1);
x=d101;
 
%figure();
%imshow(x);
  
[pop, theta] = QDEsolset_sir(npop, m, n); % Quantum-inspired Initialization
 
unfitinit = zeros(npop, 1);
for j = 1:npop
       [unfitinit(j,1)] = QDEfitnessvalue1entropy_sir(x, pop(j,:));  % Finding fitness value for each chromosome 
end
 
a=1:npop;
F1=0.25;  % Scaling Factor (DE Parameter initialization)
F2=0.25;
F3=0.2;
F4=0.2;
 
cro=0.8; % Crossover probability (DE Parameter initialization)
alpha = 0.3; % Quantum rotation angle
 
% Unified Differential Evolution Algorithm starts
qtableunfitval1(1:npop, 1:ngen) = zeros(npop, ngen);
tic; 
for i = 1:ngen
    if(i==1)
        k=1;
    else
        k=i-1;
    end
    
    [unfitsort(1:npop, i), index(:,i)]=sort(qtableunfitval1(:,k),1,'ascend');
    qtsort(1:npop,:,i)=pop(index(:,i),:,k); 
    best(:,:,i)= qtsort(end,:,i);

    for j=1:npop
        b=randperm(npop);
        c=find(b~=j);
        r1=a(b(c(1)));
        r2=a(b(c(2)));
        r3=a(b(c(3)));
        r4=a(b(c(4)));
        r5=a(b(c(5)));
        
        %[theta(j,:)] = QDErotation(theta(j,:), pop(r1,:,k), pop(r2,:,k), pop(r3,:,k), alpha);
        %dv = round(QDEmeasurement(theta(j,:), m, n));
        theta = rand() * pi;
        Q_theta = sin(theta).^2;
        dv=round(abs(pop(j,:,k)+F1*(-(best(:,:,i)-pop(j,:,k)))+ F2*(-(pop(r1,:,k)-pop(j,:,k)))+ F3*(-(pop(r2,:,k)-pop(r3,:,k)))+ F4*(-(pop(r4,:,k)-pop(r5,:,k)))+Q_theta));
        % Check and correct out-of-bound values:
        if dv(1) < 150 || dv(1) > 1350
            dv(1) = pop(r1, 1); 
        end
        if dv(2) < 150 || dv(2) > 1002
            dv(2) = pop(r1, 2);
        end
        if dv(3) < -10 || dv(3) > 10
            dv(3) = pop(r1, 3); 
        end
        if dv(4) < -10 || dv(4) > 10
            dv(4) = pop(r1, 4); 
        end
        uv = QDEopticdecross_sir(pop(j,:,k),dv,cro,m,n);
        if i > 10
            uvv = QDElocalsearch(uv, m, n, i, ngen);
        else
            uvv = uv;
        end
q1=QDEfitnessvalue1entropy_sir(x,uvv);
q2=QDEfitnessvalue1entropy_sir(x,pop(j,:,k));

if(q1 >= q2)
pop(j,:,i)= round(uvv(:,:));
qtableunfitval1(j,i)=q1;
else
pop(j,:,i) = round(pop(j,:,k));
qtableunfitval1(j,i)=q2;
end
[best_val, best_idx] = max(qtableunfitval1(:,i));
pop(1,:,i) = pop(best_idx,:,i);  
end
    
[qtableunfitval2(:,i), qtableindex(:,i)]=sort(qtableunfitval1(:,i),1,'ascend');
     
    i;
    if i == ngen
        break;
    end
end
 
t=toc;    
 
[unfitsort(1:npop, i), index(:,i)]=sort(qtableunfitval1(:,i),1,'ascend');
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
%{
figure;
plot(1:ngen, b1, '-o');
xlabel('Generation');
ylabel('Best Fitness Value');
title('Convergence Plot');
% end of part
%}
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
 
optimqtable1 = qtunfitinitsort2(end,:); 
 
groundtruth=[1180 560 10 10];
distance=sqrt(((groundtruth(1,1)-optimqtable1(1,1))^2)+((groundtruth(1,2)-optimqtable1(1,2))^2));
disp('distance:');
disp(distance);
[x1 y1 width height] = squaretable(optimqtable1);
figure();
imshow(im);
rectangle('Position',[x1,y1,width,height],'Curvature',[0.8,0.4],'LineWidth',2,'LineStyle','--');
t=toc;
fprintf('Time taken for execution: %.2f seconds\n', t);

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

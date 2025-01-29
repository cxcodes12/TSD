function [bbox, classes] = TSD_function(rgb_image)
rsz=0;
[M,N,~] = size(rgb_image);
if M*N>=800*800
    rgb_image = imresize(rgb_image,0.5);
    rsz=1;
end
 
hsv_image = rgb2hsv(rgb_image);
[M,N,~] = size(rgb_image);

% normalizare H,S,V
H = uint8(hsv_image(:,:,1) * 255);
S = uint8(hsv_image(:,:,2) * 255);
V = uint8(hsv_image(:,:,3) * 255);

% TH RED - eliminare zona acromatica, pastrare zona de rosu cromatic
mask = zeros(size(H));
for i = 1:size(H,1)
    for j = 1:size(H,2)
        if ((H(i,j) > 240 && H(i,j) <= 255) || (H(i,j) >= 0 && H(i,j) < 10))
            mask(i,j) = 255;
            if S(i,j) < 70
                 mask(i,j) = 0;
            end
%             if (V(i,j) < 50 || V(i,j) > 230)
            if (V(i,j) < 50)
                mask(i,j) = 0;
            end
        else mask(i,j) = 0;
        end
    end
end

figure(1), subplot(1,3,1), imshow(mask),title('SEG-ROSU');

% aplicare filtru median pt reducere zgomot
im_n = mask;
k1 = 5;
k2 = 5;
im_filtr = im_n;
[M,N] = size(im_n);

for i=ceil(k1/2):M-floor(k1/2)
    for j=ceil(k2/2):N-floor(k2/2)
        aux = im_n(i-floor(k1/2):i+floor(k1/2),j-floor(k2/2):j+floor(k2/2));
        aux = sort(aux(:),'ascend');
        im_filtr(i,j) = aux(ceil(length(aux)/2));
    end
end
% figure, imshow(rgb_image);
% figure, subplot(1,2,1),imshow(mask), title('imagine nefiltrata');

% aplicare operatii morfologice - inchidere (diliatare + eroziune)
img = im_filtr;
% dilatare
[M, N] = size(img);
disk = [0 1 0; 1 1 1; 0 1 0]; %origin = 1
img_b = zeros(M+2,N+2); %o img sintetica creata din img orig si o rama 0
img_dil = zeros(size(img_b));
img_b(2:M+2-1,2:N+2-1) = img;
for i=2:size(img_b,1)-1
    for j=2:size(img_b,2)-1
        patch_now = img_b(i-1:i+1, j-1:j+1);
        if(sum(sum(patch_now.*disk))~=0)
            img_dil(i,j) = 1;
        end
    end
end
img_dil = img_dil(2:end-1,2:end-1);

% eroziune
img = img_dil;
[M, N] = size(img);
disk = [0 1 0; 1 1 1; 0 1 0]; %origin = 1
img_b2 = zeros(M+2,N+2); %o img sintetica creata din img orig si o rama 0
img_er = zeros(size(img_b2));
img_b2(2:M+2-1,2:N+2-1) = img;
for i=2:M+2-1
    for j=2:N+2-1
        patch_now = img_b2(i-1:i+1, j-1:j+1);
        if (sum(sum(patch_now.*disk))==5)
            img_er(i,j) = 1;
        end
    end
end
img_er = img_er(2:end-1,2:end-1);
% subplot(1,2,2), imshow(img_er), title('img filtrata');

%% DB-SCAN pe obiectele segmentate pentru a obtine individualizarea
im = img_er;
obj = find(im==1);
col = floor((obj-1)/size(im,1))+1;
row = mod((obj-1),size(im,1))+1;
P = [col, row];

eeps = 3; %dist pt doi vectori vecini pe diag este 1.41
nPts = 10;

[cluster, noise, nr_clustere]=myDBSCAN(P,eeps,nPts);
% figure,
% for i=1:length(cluster)
%     scatter(P(cluster==i,1),P(cluster==i,2)); hold on
% end
disp(['Numarul de clustere dat de DBSCAN:',num2str(nr_clustere)]);

%% Filtrare obiecte (clustere) in functie de aspect ratio
k = 1;
clustere_ok = [];
for i=1:nr_clustere
    obj = [P(cluster==i,2),P(cluster==i,1)]; %format [linie, coloana]
    width = max(obj(:,2))-min(obj(:,2));
    height = max(obj(:,1))-min(obj(:,1));
    if (width<=2*height & width>=0.5*height & height>=round(0.05*size(mask,1)))
        clustere_ok(k) = i;
        k = k+1;
    end
end
% figure,
% for i=1:length(clustere_ok)
%     scatter(P(cluster==clustere_ok(i),1),P(cluster==clustere_ok(i),2)), hold on;
% end
% title('obiecte dupa filtrarea fct de AR');
nr_clustere_ar = length(clustere_ok);
disp(['Numarul de clustere dupa filtrarea fct de AR:',num2str(nr_clustere_ar)]);

%% pastrare obiecte (clustere) ce contin la interior si pixeli albi sau galbeni sau albastri
k=1;
clustere_ok_final = [];
for i=1:nr_clustere_ar
    c = [P(cluster==clustere_ok(i),2),P(cluster==clustere_ok(i),1)]; %format [linie, coloana]
    procent_alb = checkWhiteArea(c,S,V);
    if procent_alb>=0.1 || checkYellowArea(c,H,S,V) == 1 || checkBlueArea(c,H,S,V) == 1
        clustere_ok_final(k) = clustere_ok(i);
        k = k+1;
    end
end
% figure,
% for i=1:length(clustere_ok_final)
%     scatter(P(cluster==clustere_ok_final(i),1),P(cluster==clustere_ok_final(i),2)), hold on;
% end
% title('obiecte dupa filtrarea fct de interior alb');
nr_clustere_final = length(clustere_ok_final);
disp(['Numarul de clustere dupa filtrarea fct de interior alb:',num2str(nr_clustere_final)]);

%creare img binara cu obiectele ramase
im_obj = zeros(size(mask));
for i=1:nr_clustere_final
    contur = P(cluster==clustere_ok_final(i),:);
    for j=1:length(contur)
        im_obj(contur(j,2),contur(j,1)) = 1;   %[linie,coloana]
    end
end
figure(2), subplot(1,3,1), imshow(im_obj),title('REC-ROSU');

% figure, imshow(im_obj);
% imwrite(im_obj,'test99_maskred_eliminare.png');

%% Analiza fiecarui ROI
clc
classes = {};
bbox = [];
for k=1:nr_clustere_final
    c = [P(cluster==clustere_ok_final(k),2), P(cluster==clustere_ok_final(k),1)]; %format [linie, coloana]
    im_c = zeros(size(mask));
    for i=1:size(c,1)
        im_c(c(i,1),c(i,2)) = 1;
    end
    mask2 = im_c;

    %obtinere coordonate BB ROI
    im_c = mask2;
    [rows,cols] = find(im_c==1);
    c = [rows,cols];
    row = [min(c(:,1)), max(c(:,1))];
    col = [min(c(:,2)), max(c(:,2))];
    x = col(1);
    y = row(1);
    width = col(2)-col(1);
    height = row(2) - row(1);
    
    %crearea unei margini negre in jurul ROI
    im_c = im_c(row(1):row(2),col(1):col(2));
    bs = 7;
    im_b = zeros(size(im_c,1)+2*bs, size(im_c,2)+2*bs);
    im_b(bs+1:end-bs,bs+1:end-bs) = im_c;
    im_c = im_b; 
    
    %umplere goluri in img binara ROI
    for i=1:size(im_c,1)
        line = im_c(i,:);
        if ~isempty(find(line==1))
            idx = [min(find(line==1)), max(find(line==1))];
            im_c(i,idx(1):idx(2)) = 1;
        end
    end
    
    
    %verificarea conditiilor pentru SQUARE/TDP/TUP
    [M,N] = size(im_c);
    if M*N <=60*60
        im_c = imresize(im_c,[2*M,2*N]);
    end
    %aplicare filtru gaussian
    msigma = 0.5;
    im_c_f = myGaussianFilter(im_c,msigma);
    %aplicare Harris Corner Detector pe cele 6 zone prestabilite ale imaginii
    [check, shape] = checkCorners(im_c_f);
    if check==1 && (strcmp(shape,'Cedeaza')||strcmp(shape,'Triunghi'));
        bbox = [bbox; x,y,width,height];
        if strcmp(shape,'Cedeaza')
            classes = vertcat(classes,'Cedeaza');
        else classes = vertcat(classes,'Avertizare');
        end
    else
    %verificare simetrie radiala (FRS)

    % parametri pentru algoritmul FRS
    r_col = floor((col(2)-col(1))/2);
    r_row = floor((row(2)-row(1))/2);

    raza = max(r_col,r_row); % raza pentru calculul simetriei radiale
    raza = raza - round(0.1*raza); 
    nr_unghiuri = 60; % nr de unghiuri luate in calcul
    TH = 1; % fiind o imagine binara, TH este 1

    [frs_decision,outputimage] = myFRS(im_c, raza, nr_unghiuri, TH);
    [Ex,Oct] = myEOmeasures(im_c);
    if frs_decision==1 && Ex > Oct && Ex > 0.9
        shape = 'circle';
        classes = vertcat(classes,'Interzicere');
        bbox = [bbox; x,y,width,height];
    else if frs_decision==1 && Ex < Oct && Oct > 0.9
            shape = 'octagon';
            classes = vertcat(classes, 'STOP');
            bbox = [bbox; x,y,width,height];
        end
    end
    end
end


%% TH BLUE - eliminare zona acromatica, pastrare zona de albastru cromatic
mask = zeros(size(H));
for i = 1:size(H,1)
    for j = 1:size(H,2)
        if (H(i,j) >= 130 && H(i,j) <= 180)
            mask(i,j) = 255;
            if S(i,j) < 130
                 mask(i,j) = 0;
            end
            if (V(i,j) < 65)
                mask(i,j) = 0;
            end
        else mask(i,j) = 0;
        end
    end
end
figure(1), subplot(1,3,2), imshow(mask),title('ALBASTRU');

% aplicare filtru median pt reducere zgomot
im_n = mask;
k1 = 5;
k2 = 5;
im_filtr = im_n;
[M,N] = size(im_n);

for i=ceil(k1/2):M-floor(k1/2)
    for j=ceil(k2/2):N-floor(k2/2)
        aux = im_n(i-floor(k1/2):i+floor(k1/2),j-floor(k2/2):j+floor(k2/2));
        aux = sort(aux(:),'ascend');
        im_filtr(i,j) = aux(ceil(length(aux)/2));
    end
end
% figure, imshow(rgb_image);
% figure, subplot(1,2,1),imshow(mask), title('imagine nefiltrata');

% aplicare operatii morfologice - inchidere (diliatare + eroziune)
img = im_filtr;
% dilatare
[M, N] = size(img);
disk = [0 1 0; 1 1 1; 0 1 0]; %origin = 1
img_b = zeros(M+2,N+2); %o img sintetica creata din img orig si o rama 0
img_dil = zeros(size(img_b));
img_b(2:M+2-1,2:N+2-1) = img;
for i=2:size(img_b,1)-1
    for j=2:size(img_b,2)-1
        patch_now = img_b(i-1:i+1, j-1:j+1);
        if(sum(sum(patch_now.*disk))~=0)
            img_dil(i,j) = 1;
        end
    end
end
img_dil = img_dil(2:end-1,2:end-1);

% eroziune
img = img_dil;
[M, N] = size(img);
disk = [0 1 0; 1 1 1; 0 1 0]; %origin = 1
img_b2 = zeros(M+2,N+2); %o img sintetica creata din img orig si o rama 0
img_er = zeros(size(img_b2));
img_b2(2:M+2-1,2:N+2-1) = img;
for i=2:M+2-1
    for j=2:N+2-1
        patch_now = img_b2(i-1:i+1, j-1:j+1);
        if (sum(sum(patch_now.*disk))==5)
            img_er(i,j) = 1;
        end
    end
end
img_er = img_er(2:end-1,2:end-1);
% subplot(1,2,2), imshow(img_er), title('img filtrata');

%% DB-SCAN pe obiectele segmentate pentru a obtine individualizarea
im = img_er;
obj = find(im==1);
col = floor((obj-1)/size(im,1))+1;
row = mod((obj-1),size(im,1))+1;
P = [col, row];

eeps = 3.5; %dist pt doi vectori vecini pe diag este 1.41
nPts = 10;

[cluster, noise, nr_clustere]=myDBSCAN(P,eeps,nPts);
% figure,
% for i=1:length(cluster)
%     scatter(P(cluster==i,1),P(cluster==i,2)); hold on
% end
disp(['Numarul de clustere dat de DBSCAN:',num2str(nr_clustere)]);

%% Filtrare obiecte (clustere) in functie de aspect ratio
k = 1;
clustere_ok = [];
for i=1:nr_clustere
    obj = [P(cluster==i,2),P(cluster==i,1)]; %format [linie, coloana]
    width = max(obj(:,2))-min(obj(:,2));
    height = max(obj(:,1))-min(obj(:,1));
    if (width<=5*height & width>=0.5*height & height>=round(0.05*size(mask,1)))
        clustere_ok(k) = i;
        k = k+1;
    end
end
% figure,
% for i=1:length(clustere_ok)
%     scatter(P(cluster==clustere_ok(i),1),P(cluster==clustere_ok(i),2)), hold on;
% end
% title('obiecte dupa filtrarea fct de AR');
nr_clustere_ar = length(clustere_ok);
disp(['Numarul de clustere dupa filtrarea fct de AR:',num2str(nr_clustere_ar)]);

%% pastrare obiecte (clustere) ce contin la interior si pixeli albi
k=1;
clustere_ok_final = [];
for i=1:nr_clustere_ar
    c = [P(cluster==clustere_ok(i),2),P(cluster==clustere_ok(i),1)]; %format [linie, coloana]
    if checkWhiteArea(c,S,V)>=0.1
        clustere_ok_final(k) = clustere_ok(i);
        k = k+1;
    end
end
% figure,
% for i=1:length(clustere_ok_final)
%     scatter(P(cluster==clustere_ok_final(i),1),P(cluster==clustere_ok_final(i),2)), hold on;
% end
% title('obiecte dupa filtrarea fct de interior alb');
nr_clustere_final = length(clustere_ok_final);
%creare img binara cu obiectele ramase
im_obj = zeros(size(mask));
for i=1:nr_clustere_final
    contur = P(cluster==clustere_ok_final(i),:);
    for j=1:length(contur)
        im_obj(contur(j,2),contur(j,1)) = 1;   %[linie,coloana]
    end
end
figure(2), subplot(1,3,2), imshow(im_obj),title('REC-ALBASTRU');

nr_clustere_final = length(clustere_ok_final);
disp(['Numarul de clustere dupa filtrarea fct de interior alb:',num2str(nr_clustere_final)]);

%% Analiza fiecarui ROI
clc
for k=1:nr_clustere_final
    c = [P(cluster==clustere_ok_final(k),2), P(cluster==clustere_ok_final(k),1)]; %format [linie, coloana]
    im_c = zeros(size(mask));
    for i=1:size(c,1)
        im_c(c(i,1),c(i,2)) = 1;
    end
    mask2 = im_c;

    %obtinere coordonate BB ROI
    im_c = mask2;
    [rows,cols] = find(im_c==1);
    c = [rows,cols];
    row = [min(c(:,1)), max(c(:,1))];
    col = [min(c(:,2)), max(c(:,2))];
    x = col(1);
    y = row(1);
    width = col(2)-col(1);
    height = row(2) - row(1);
    
    %crearea unei margini negre in jurul ROI
    im_c = im_c(row(1):row(2),col(1):col(2));
    bs = 5;
    im_b = zeros(size(im_c,1)+2*bs, size(im_c,2)+2*bs);
    im_b(bs+1:end-bs,bs+1:end-bs) = im_c;
    im_c = im_b; 
    
    %umplere goluri in img binara ROI
    for i=1:size(im_c,1)
        line = im_c(i,:);
        if ~isempty(find(line==1))
            idx = [min(find(line==1)), max(find(line==1))];
            im_c(i,idx(1):idx(2)) = 1;
        end
    end

    
    %verificarea conditiilor pentru SQUARE/TDP/TUP
    
    %aplicare filtru gaussian
    [M,N] = size(im_c);
    if M*N <=60*60
        im_c = imresize(im_c,[2*M,2*N]);
    end
    msigma = 0.3;
    im_c_f = myGaussianFilter(im_c,msigma);
%   aplicare Harris Corner Detector pe cele 6 zone prestabilite ale imaginii
    [check, shape] = checkCorners(im_c_f);
    if check==1 && (strcmp(shape,'Patrat'))
        bbox = [bbox; x,y,width,height];
        classes = vertcat(classes,'Orientare_Info');
    else
    %verificare simetrie radiala (FRS) si Excentricitate / Octogonalitate

    % parametri pentru algoritmul FRS
    r_col = floor((col(2)-col(1))/2);
    r_row = floor((row(2)-row(1))/2);

    raza = max(r_col,r_row); % raza pentru calculul simetriei radiale
    raza = raza - round(0.1*raza); 
    nr_unghiuri = 60; % nr de unghiuri luate in calcul
    TH = 1; % fiind o imagine binara, TH este 1

    [frs_decision,outputimage] = myFRS(im_c, raza, nr_unghiuri, TH);
    [Ex] = myEOmeasures(im_c);
    if frs_decision==1 && Ex > 0.9
        shape = 'circle';
        classes = vertcat(classes, 'Obligare');
        bbox = [bbox; x,y,width,height];
    end
    end
end

%% TH YELLOW - eliminare zona acromatica, pastrare zona de galben cromatic - prioritate

mask = zeros(size(H));
for i = 1:size(H,1)
    for j = 1:size(H,2)
        if (H(i,j) >= 15 && H(i,j) <= 45)
            mask(i,j) = 255;
            if S(i,j) < 70
                 mask(i,j) = 0;
            end
            if (V(i,j) < 50)
                mask(i,j) = 0;
            end
        else mask(i,j) = 0;
        end
    end
end
figure(1), subplot(1,3,3), imshow(mask),title('SEG-GALBEN');

% aplicare filtru median pt reducere zgomot
im_n = mask;
k1 = 5;
k2 = 5;
im_filtr = im_n;
[M,N] = size(im_n);

for i=ceil(k1/2):M-floor(k1/2)
    for j=ceil(k2/2):N-floor(k2/2)
        aux = im_n(i-floor(k1/2):i+floor(k1/2),j-floor(k2/2):j+floor(k2/2));
        aux = sort(aux(:),'ascend');
        im_filtr(i,j) = aux(ceil(length(aux)/2));
    end
end
% figure, imshow(rgb_image);
% figure, subplot(1,2,1),imshow(mask), title('imagine nefiltrata');

% aplicare operatii morfologice - inchidere (diliatare + eroziune)
img = im_filtr;
% dilatare
[M, N] = size(img);
disk = [0 1 0; 1 1 1; 0 1 0]; %origin = 1
img_b = zeros(M+2,N+2); %o img sintetica creata din img orig si o rama 0
img_dil = zeros(size(img_b));
img_b(2:M+2-1,2:N+2-1) = img;
for i=2:size(img_b,1)-1
    for j=2:size(img_b,2)-1
        patch_now = img_b(i-1:i+1, j-1:j+1);
        if(sum(sum(patch_now.*disk))~=0)
            img_dil(i,j) = 1;
        end
    end
end
img_dil = img_dil(2:end-1,2:end-1);

% eroziune
img = img_dil;
[M, N] = size(img);
disk = [0 1 0; 1 1 1; 0 1 0]; %origin = 1
img_b2 = zeros(M+2,N+2); %o img sintetica creata din img orig si o rama 0
img_er = zeros(size(img_b2));
img_b2(2:M+2-1,2:N+2-1) = img;
for i=2:M+2-1
    for j=2:N+2-1
        patch_now = img_b2(i-1:i+1, j-1:j+1);
        if (sum(sum(patch_now.*disk))==5)
            img_er(i,j) = 1;
        end
    end
end
img_er = img_er(2:end-1,2:end-1);
% subplot(1,2,2), imshow(img_er), title('img filtrata');

%% DB-SCAN pe obiectele segmentate pentru a obtine individualizarea
im = img_er;
obj = find(im==1);
col = floor((obj-1)/size(im,1))+1;
row = mod((obj-1),size(im,1))+1;
P = [col, row];

eeps = 2; %dist pt doi vectori vecini pe diag este 1.41
nPts = 12;

[cluster, noise, nr_clustere]=myDBSCAN(P,eeps,nPts);
% figure,
% for i=1:length(cluster)
%     scatter(P(cluster==i,1),P(cluster==i,2)); hold on
% end
disp(['Numarul de clustere dat de DBSCAN:',num2str(nr_clustere)]);

%% Filtrare obiecte (clustere) in functie de aspect ratio
k = 1;
clustere_ok = [];
for i=1:nr_clustere
    obj = [P(cluster==i,2),P(cluster==i,1)]; %format [linie, coloana]
    width = max(obj(:,2))-min(obj(:,2));
    height = max(obj(:,1))-min(obj(:,1));
    if (width<=2*height & width>=0.5*height  & height>=round(0.05*size(mask,1)))
        clustere_ok(k) = i;
        k = k+1;
    end
end
% figure,
% for i=1:length(clustere_ok)
%     scatter(P(cluster==clustere_ok(i),1),P(cluster==clustere_ok(i),2)), hold on;
% end
% title('obiecte dupa filtrarea fct de AR');
nr_clustere_ar = length(clustere_ok);
disp(['Numarul de clustere dupa filtrarea fct de AR:',num2str(nr_clustere_ar)]);

%% pastrare obiecte (clustere) in functie de omogenitate
k=1;
clustere_ok_final = [];
for t=1:nr_clustere_ar
    checkH = [];
    c = [P(cluster==clustere_ok(t),2),P(cluster==clustere_ok(t),1)]; %format [linie, coloana]
    im_ob = zeros(size(S));
    for i=1:size(c,1)
      im_ob(c(i,1),c(i,2)) = 1;
    end
    for i=1:size(im_ob,1)
        line = im_ob(i,:);
        if ~isempty(find(line==1))
            col = [min(find(line==1)), max(find(line==1))];
            for j=col(1):col(2)
                checkH = [checkH, H(i,j)];
            end
        end
    end
    if sum(checkH>=15 & checkH<=45)/numel(checkH)>=0.95
        clustere_ok_final(k) = clustere_ok(t);
        k = k+1;
    end
end


% figure,
% for i=1:length(clustere_ok_final)
%     scatter(P(cluster==clustere_ok_final(i),1),P(cluster==clustere_ok_final(i),2)), hold on;
% end
% title('obiecte dupa filtrarea fct de omogenitate');
nr_clustere_final = length(clustere_ok_final);

im_obj = zeros(size(mask));
for i=1:nr_clustere_final
    contur = P(cluster==clustere_ok_final(i),:);
    for j=1:length(contur)
        im_obj(contur(j,2),contur(j,1)) = 1;   %[linie,coloana]
    end
end
figure(2), subplot(1,3,3), imshow(im_obj),title('REC-GALBEN');

% figure, imshow(im_obj);
% imwrite(im_obj,'test99_maskyellow_eliminare.png');

%% Analiza fiecarui ROI
clc
for k=1:nr_clustere_final
    c = [P(cluster==clustere_ok_final(k),2), P(cluster==clustere_ok_final(k),1)]; %format [linie, coloana]
    im_c = zeros(size(mask));
    for i=1:size(c,1)
        im_c(c(i,1),c(i,2)) = 1;
    end
    mask2 = im_c;
    
    %obtinere coordonate BB ROI
    im_c = mask2;
    [rows,cols] = find(im_c==1);
    row = [min(c(:,1)), max(c(:,1))];
    col = [min(c(:,2)), max(c(:,2))];
    x = col(1);
    y = row(1);
    width = col(2)-col(1);
    height = row(2) - row(1);
    
    %crearea unei margini negre in jurul ROI
    im_c = im_c(row(1):row(2),col(1):col(2));
    bs = 5;
    im_b = zeros(size(im_c,1)+2*bs, size(im_c,2)+2*bs);
    im_b(bs+1:end-bs,bs+1:end-bs) = im_c;
    im_c = im_b; 
    
    %umplere goluri in img binara ROI
    for i=1:size(im_c,1)
        line = im_c(i,:);
        if ~isempty(find(line==1))
            idx = [min(find(line==1)), max(find(line==1))];
            im_c(i,idx(1):idx(2)) = 1;
        end
    end
    
    im_c = imresize(im_c,[2*size(im_c,1),2*size(im_c,2)]);
    
    % sa adaug aici eroziune sau dilatare, ceva?
    
    %verificarea conditiilor pentru SQUARE/TDP/TUP
    
    %aplicare filtru gaussian
    im_c = imrotate(im_c,45);
    msigma = 0.3;
    im_c_f = myGaussianFilter(im_c,msigma);
    %aplicare Harris Corner Detector pe cele 6 zone prestabilite ale imaginii
    [row,col] = find(im_c_f~=0);
    im_c_f = im_c_f(min(row)-3:max(row)+3,min(col)-3:max(col)+3);
    [check, shape] = checkCorners(im_c_f);
    if check==1 && (strcmp(shape,'Patrat'))
        bbox = [bbox; x,y,width,height];
        classes = vertcat(classes,'Prioritate');
    else
        
    end
end
if rsz
bbox = 2*bbox;
end
end

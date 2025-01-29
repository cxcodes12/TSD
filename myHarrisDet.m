function [isCorner, rez] = myHarrisDet(im)
isCorner = 0;
isFlat = 0;
isEdge = 0;
dx = [-1, 0, 1];
dy = [-1; 0; 1];

%calcul gradienti pe x si pe y
% im = imresize(im,[2*size(im,1),2*size(im,2)]);
[M,N] = size(im);
Ix = zeros(M,N);
for i=1:M
    for j=2:N-1
        Ix(i,j) = sum(sum(dx.*im(i,j-1:j+1)));
    end
end
Iy = zeros(M,N);
for i=2:M-1
    for j=1:N
        Iy(i,j) = sum(sum(dy.*im(i-1:i+1,j)));
    end
end
Ixy = Ix.*Iy;

%parcurg cu o fereastra 2k+1 x 2k+1 si calculez parametrii matricei M pt fiecare pct 
%Matricea M - matricea de autocorelatie (matricea de moment 2)
R = zeros(M,N);
k = 5;
for i=k+1:M-k
    for j=k+1:N-k
        Ix2 = sum(sum(Ix(i-k:i+k,j-k:j+k).^2));
        Iy2 = sum(sum(Iy(i-k:i+k,j-k:j+k).^2));
        Ixy2 = sum(sum(Ixy(i-k:i+k,j-k:j+k)));
        M = [Ix2, Ixy2; Ixy2, Iy2];
        kv = 0.05; %definit empiric in intervalul 0.04 - 0.06
        R(i,j) = det(M) - kv*(trace(M))^2;
    end
end

if max(R(:))<0
    isEdge = 1; 
end
if max(R(:)) < 50 & R>=0
    isFlat = 1;
end
if max(R(:)) > 50
    isCorner = 1;
end

rez = max(R(:));
% de modificat aceste praguri in caz de erori majore

end

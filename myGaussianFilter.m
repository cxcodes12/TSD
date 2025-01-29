function im_filtr = myGaussianFilter(im, msigma)
%msigma = abaterea medie, la patrat devine dispersie, cu cat dispersia e mai mare, cu atat banda e mai ingusta

%creare suprafata gaussiana (masca)
x = 0:0.5:1;
y = 0:0.5:1;
G = zeros(length(x),length(y));
[k1,k2] = size(G);

%adaugare rama astfel incat sa am imaginea neafectata de efectul convolutiei

im_b = zeros(size(im,1)+2*k1,size(im,2)+2*k2);
im_b((k1+1):end-k1,(k2+1):end-k2) = im;
im = im_b;

%calcul gaussiana
for i=1:length(x)
    for j=1:length(y)
        G(i,j) = 1/(2*pi*msigma^2)*exp(-(x(i)^2+y(j)^2)/(2*msigma^2));
    end
end
G = G/sum(G(:));  %normare

%aplicare masca prin convolutie
[k1,k2] = size(G);
[M,N] = size(im);
im_filtr = zeros(M,N);
for i=ceil(k1/2):M-floor(k1/2)
    for j=ceil(k2/2):N-floor(k2/2)
        im_filtr(i,j) = sum(sum(double(im(i-floor(k1/2):i+floor(k1/2),j-floor(k2/2):j+floor(k2/2))).*G));
    end
end
im_filtr = im_filtr((k1+1):M-k1,(k2+1):N-k1);

end

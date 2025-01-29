function [E,O] = myEOmeasures(im)

[M,N] = size(im);
[x, y] = meshgrid(1:N,1:M);

M01 = sum(sum(y.*im));
M10 = sum(sum(x.*im));
M00 = sum(sum(im));
centroid_X = M10 / M00;
centroid_Y = M01 / M00;

mu00 = sum(im(:));
mu11 = sum((x(:) - centroid_X).*(y(:) - centroid_Y).*im(:));
mu02 = sum((y(:)-centroid_Y).^2.*im(:));
mu20 = sum((x(:)-centroid_X).^2.*im(:));

I1 = (mu20*mu02-mu11^2)/(mu00^4);
if I1 <= 1/(16*(pi^2))
    E = 16*(pi^2)*I1;
else E = 1/(16*(pi^2)*I1);
end

if I1 <= 1/(15.932*(pi^2))
    O = 15.932*(pi^2)*I1;
else O = 1/(15.932*(pi^2)*I1);
end

end

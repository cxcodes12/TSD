function [blueArea, nr_albastru] = checkBlueArea(c,H,S,V)
%pentru fiecare linie din forma obiectului c, aflu minimul si maximul pe
%coloane, iar intre cele doua coloane limite, numar cati pixeli de
%albastru sunt
im_ob = zeros(size(S));
for i=1:size(c,1)
    im_ob(c(i,1),c(i,2)) = 1;
end

nr_albastru = 0; %nr de pixeli de albastru
area = 0; %nr total de pixeli ai obiectului incluzand interiorul
for i=1:size(im_ob,1)
    line = im_ob(i,:);
    if ~isempty(find(line==1))
        col = [min(find(line==1)), max(find(line==1))];
        for j=col(1):col(2)
            if(H(i,j)>=130 && H(i,j)<=180 && S(i,j)>=130 && V(i,j)>= 65)
                nr_albastru= nr_albastru+1;
            end
            area = area+1;
        end
    end
end
if(nr_albastru / area >= 0.2)
    blueArea = 1;       
else blueArea = 0;


end
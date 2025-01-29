function [checkYellow, nr_galben] = checkYellowArea(c,H,S,V)
%pentru fiecare linie din forma obiectului c, aflu minimul si maximul pe
%coloane, iar intre cele doua coloane limite, numar cati pixeli de
%galben sunt
im_ob = zeros(size(S));
for i=1:size(c,1)
    im_ob(c(i,1),c(i,2)) = 1;
end

nr_galben = 0; %nr de pixeli de galben
area = 0; %nr total de pixeli ai obiectului incluzand interiorul
for i=1:size(im_ob,1)
    line = im_ob(i,:);
    if ~isempty(find(line==1))
        col = [min(find(line==1)), max(find(line==1))];
        for j=col(1):col(2)
            if(H(i,j)>=20 && H(i,j)<=40 && S(i,j)>=160 && V(i,j)>= 100)
                nr_galben = nr_galben+1;
            end
            area = area+1;
        end
    end
end
if(nr_galben / area >= 0.2)
    checkYellow = 1;       
else checkYellow = 0;


end
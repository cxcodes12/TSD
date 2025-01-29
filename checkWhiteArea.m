function [procent_alb] = checkWhiteArea(c,S,V)
%pentru fiecare linie din forma obiectului c, aflu minimul si maximul pe
%coloane, iar intre cele doua coloane limite, numar cati pixeli de
%alb sunt
im_ob = zeros(size(S));
for i=1:size(c,1)
    im_ob(c(i,1),c(i,2)) = 1;
end

nr_alb = 0; %nr de pixeli de alb
area = 0; %nr total de pixeli ai obiectului incluzand interiorul
for i=1:size(im_ob,1)
    line = im_ob(i,:);
    if ~isempty(find(line==1))
        col = [min(find(line==1)), max(find(line==1))];
        for j=col(1):col(2)
            if(S(i,j)<120 && V(i,j)>80)
                nr_alb = nr_alb+1;
            end
            area = area+1;
        end
    end
end
procent_alb = nr_alb/area;



end
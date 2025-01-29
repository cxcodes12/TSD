function [check, shape] = checkCorners(im_c)

%initiere parametri de masurare a probabilitatii
TL=0; TC=0; TR=0; BL=0; BC=0; BR=0;
[M,N] = size(im_c);
% Top Left
im = im_c(1:round(M/4),1:round(N/4));
[isCorner_TL, RTL] = myHarrisDet(im);
if isCorner_TL
    TL = 0.25;
end
% Top Centre
im = im_c(1:round(M/4),round(N/2-N/8):round(N/2+N/8));
[isCorner_TC, RTC] = myHarrisDet(im);
if isCorner_TC
    TC = 0.34;
end
% Top Right
im = im_c(1:round(M/4),round(3*N/4):N);
[isCorner_TR, RTR] = myHarrisDet(im);
if isCorner_TR
    TR = 0.25;
end
% Bottom Left
im = im_c(round(3*M/4):M,1:round(N/4));
[isCorner_BL, RBL] = myHarrisDet(im);
if isCorner_BL
    BL = 0.25;
end
% Bottom Centre
im = im_c(round(3*M/4):M,round(N/2-N/8):round(N/2+N/8));
[isCorner_BC, RBC] = myHarrisDet(im);
if isCorner_BC
    BC = 0.34;
end
% Bottom Right
im = im_c(round(3*M/4):M,round(3*N/4):N);
[isCorner_BR, RBR] = myHarrisDet(im);
if isCorner_BR
    BR = 0.25;
end

% probabilitati pentru formele:
sqp = TL+TR+BL+BR; %patrat sau dreptunghi
tdp = 1.32*(BL+BR)+TC-1.1*(TL+TR); %triunghi atentionare
tup = 1.32*(TL+TR)+BC-1.1*(BL+BR); %triunghi cedeaza
shapes_p = [sqp,tdp,tup];
rez = max(shapes_p);
check = 0;
if rez>=0.6
    if rez==sqp
        shape = 'Patrat';
        check = 1;
    end
    if rez==tdp
        shape = 'Triunghi';
        check = 1;
    end
    if rez==tup
        shape = 'Cedeaza';
        check = 1;
    end
else check = 0; shape = 'none';
end


end
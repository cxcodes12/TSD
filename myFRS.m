function [checkCircle,output_image] = myFRS(input_image, r, num_angles, th)
    output_image = zeros(size(input_image));
    %parcurgere pixeli si verificare diferenta pixeli la distanta radius
    for i = 1:size(input_image, 1)
        for j = 1:size(input_image, 2)
            symmetries = zeros(1, num_angles);

            for k = 1:num_angles
                angle = pi * k / num_angles;
                x1 = round(j + r * cos(angle));
                y1 = round(i + r * sin(angle));
                x2 = round(j - r * cos(angle));
                y2 = round(i - r * sin(angle));

                if x1 >= 1 && x1 <= size(input_image, 2) && y1 >= 1 && y1 <= size(input_image, 1) && x2 >= 1 && x2 <= size(input_image, 2) && y2 >= 1 && y2 <= size(input_image, 1)
                    symmetries(k) = input_image(y1,x1) + input_image(y2,x2);
                end
            end

            output_image(i,j) = round(255*sum(symmetries==2*th)/num_angles);
        end
    end
    %verificarea maximului din centrul imaginii
    [M,N] = size(output_image);
    centre = [round(M/2), round(N/2)];
    p = round(0.1*max(M,N));
    im = output_image(centre(1)-p:centre(1)+p, centre(2)-p:centre(2)+p);
    if max(im(:)) == max(output_image(:)) & max(im(:))>=150
        checkCircle = 1;
    else checkCircle = 0;
    
    
end
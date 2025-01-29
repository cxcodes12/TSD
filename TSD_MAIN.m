%main code for traffic signs detection in an image
clc
clearvars
close all

im = imread('type_your_image_name');
[bboxes,classes] = TSD_function(im);

if (size(bboxes,1)~=0)
    figure, imshow(im), hold on
    for i = 1:size(bboxes,1)
        rectangle('Position', [bboxes(i,1),bboxes(i,2), bboxes(i,3), bboxes(i,4)], 'EdgeColor', 'g', 'LineWidth', 4); 
        text(bboxes(i,1), bboxes(i,2) - 10, classes(i), 'Color', 'm', 'FontSize', 18, 'FontWeight', 'bold'); 
        hold on
    end
else disp('no traffic signs detected');
end





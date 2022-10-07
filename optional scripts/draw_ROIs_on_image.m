%% Create ROIs on widefield image
% Draw ROIs on widefield image. Select as many ROIs as
% needed, by clicking draw (or enter) after each ROI. Click on quit at the
% end. Possible to zoom on the image, to create more precise ROIs.
% 
% Output is similar to PSD mask script, cell with all the boundary
% coordinates of the ROIs in nm.
% 
% Manon Westra, April 2021

clear; clc; close all;
draw_ROIs_on_image_param            % load parameters

% load/display image
imagelist = dir(filename);
widefield = imread(imagelist.name);

fig = figure('name','Widefield to draw ROIs on');

% correct origin image and axes image to nm scale
RI = imref2d(size(widefield));
RI.XWorldLimits = [(0-(pixel/2)) (size(widefield,2)*pixel-(pixel/2))];
RI.YWorldLimits = [(0-(pixel/2)) (size(widefield,1)*pixel-(pixel/2))];

% plot widefield with correct axes
imshow(widefield,RI,[low high]);axis image;

% select ROIs for spines
again = true;
i = 0;
while again
    promptMessage = sprintf('Draw SPINE ROI #%d or Quit?', i + 1);
	titleBarCaption = 'Continue?';
	button = questdlg(promptMessage, titleBarCaption, 'Draw', 'Quit', 'Draw');
    if strcmpi(button, 'Quit')
		break;
    end
    i = i + 1;
    ROIoutline = drawfreehand;
    spines{i,:}(:,2) = ROIoutline.Position(:,1); % x coordinates similar as psdmask in the second column
    spines{i,:}(:,1) = ROIoutline.Position(:,2); % y coordinates similar as psdmask in the first column
    
end

disp('ROIs for spines created!');

% select ROIs for dendrites
again = true;
i = 0;
while again
    promptMessage = sprintf('Draw DENDRITE ROI #%d or Quit?', i + 1);
	titleBarCaption = 'Continue?';
	button = questdlg(promptMessage, titleBarCaption, 'Draw', 'Quit', 'Draw');
    if strcmpi(button, 'Quit')
		break;
    end
    i = i + 1;
    ROIoutline = drawfreehand;
    dendrites{i,:}(:,2) = ROIoutline.Position(:,1); % x coordinates similar as psdmask in the second column
    dendrites{i,:}(:,1) = ROIoutline.Position(:,2); % y coordinates similar as psdmask in the first column
    
end

% saving
saveas(gcf,'SpinesDendrites');
save spinedendrite.mat spines dendrites
disp('ROIs for dendrites created and data saved!');
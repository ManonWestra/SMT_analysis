function h = PlotBoundaries(wf_image, pixel)
% plot outline widefield image using mask
% Manon Westra, 2022

Igray = imread(wf_image);
Ibw = imbinarize(Igray,'adaptive');
Ifill = imfill(Ibw,'holes');
[B,L] = bwboundaries(Ifill);
B = cellfun(@(x) x*pixel,B,'un',0);
visboundaries(B, 'Color', 'k', 'LineStyle', ':','LineWidth', 0.5)

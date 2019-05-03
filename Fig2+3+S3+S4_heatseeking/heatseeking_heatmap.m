function peak = heatseeking_heatmap(pos,path,name,ifdot,cmax)

% 12x16 bins
ctrs = cell(2,1);
ctrs{1} = 20:40:640;
ctrs{2} = 20:40:480;

% test plot of the bins
% figure,hold on;
% imshow(dBR);
% for i=1:length(ctrs{1})
%     text(ctrs{1}(i),240,num2str(i),'Color','r');
% end
% for j=1:length(ctrs{2})
%     text(320,ctrs{2}(j),num2str(j),'Color','b');
% end
% hold off;

% same orientation as the picture
H = rot90(hist3(pos,ctrs));
%H = H / sum(sum(H));
peak = max(max(H));

%plots with Peltier array
% Peltier = [472,125;475,421;135,404;134,149];
Peltier = [1.369011799410029e+02,1.512964601769910e+02;...
    4.682286135693216e+02,1.352492625368730e+02;...
    4.710604719764012e+02,4.118274336283185e+02;...
    1.425648967551622e+02,3.957802359882005e+02];
Peltier(:,2) = 480 - Peltier(:,2);
Peltier = (Peltier ./ 40) + 0.5; % scales to the bins
h = figure; hold on;
contourf(H,16,'LineStyle','none'),shading flat;
colormap(flipud(gray));
if nargin == 5
    caxis([0 cmax]);
end
impoly(gca,Peltier);
if ifdot
    load([path 'ellipse.mat'],'polypos');
    polypos(:,2) = 480 - polypos(:,2);
    polypos = (polypos ./ 40) + 0.5;
    impoly(gca,polypos);
end
set(gca,'DataAspectRatio',[1 1 1]);
xlim([0.5,16.5]);
ylim([0.5,12.5]);
hold off;

if ~isempty(path) && ~isempty(name)
    saveas(h,[path '/' name '.pdf']);
end

end
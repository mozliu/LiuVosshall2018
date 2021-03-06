path = 'D:\Dropbox\labwork\behavior\heatseeking-pilot\';
dot = {'20171011-orl1003-1',...
    '20171014-orl1003-0',...
    '20171016-orl1003-1',...
    '20171019-orl1010-0',...
    '20171019-orl1010-2',...
      '20171021_orl1010-1',...
      '20171021_orl1010-3',...
      '20171023_orl1010-1'};
blank = {'20171014-orl1003-1',...
    '20171016-orl1003-0',...
    '20171016-orl1003-2',...
    '20171019-orl1010-1',...
    '20171021_orl1010-0',...
        '20171021_orl1010-2',...
        '20171023_orl1010-0',...
        '20171023_orl1010-2'};
  
templabels = [400,33.5,26,55,45,60,31,38.5,40,28.5,36,50];
dotcmax = zeros(1,12);

[allpos,cdot,cdotn] = gather_all_pos(blank,path,'shuffle_ellipse.mat');
for i=1:12
    %dotcmax(i) = heatseeking_heatmap(allpos{i},'','',1);
    heatseeking_heatmap(allpos{i},path,['hm_shufblank_' num2str(templabels(i)) '.pdf'],0,360);
    title(['orl-' num2str(templabels(i)) '.pdf']);
    close
end
[allpos,ddot,ddotn] = gather_all_pos(dot,path,'shuffle_ellipse.mat');
for i=1:12
    %dotcmax(i) = heatseeking_heatmap(allpos{i},'','',1);
    heatseeking_heatmap(allpos{i},path,['hm_shufdot_' num2str(templabels(i)) '.pdf'],1,360);
    title(['dot-' num2str(templabels(i)) '.pdf']);
    %close
end

save(fullfile(path,'shuffle_dot'),'cdot','cdotn','ddot','ddotn');
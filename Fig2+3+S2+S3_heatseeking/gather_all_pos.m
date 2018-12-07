function [allpos,all_on_dot,all_dot_n] = gather_all_pos(fn,path,ellipse_name)

if nargin < 3
    ellipse_name = 'ellipse.mat';
    if nargin < 2
        path = 'D:\Dropbox\labwork\behavior\heatseeking-pilot\';
    end
end
% orl = {'20180627_2',...
%        '20180628_2',...
%        '20180629_3',...
%        '20180711_1',...
%        '20180712_1',...
%        '20180713_2'};
% orlgr3 = {'20180627_3',...
%        '20180628_1',...
%        '20180629_1',...
%        '20180711_2',...
%        '20180712_2',...
%        '20180713_1'};
% gr3 = {'20180627_1',...
%        '20180628_3',...
%        '20180629_2',...
%        '20180711_3',...
%        '20180712_3',...
%        '20180713_3'};

% declare coordinates to draw on

load(fullfile(path,ellipse_name),'dotpos','polypos');
rng('shuffle');

% when to sample
interval = 90;
sampling = 630:720:9180;

allpos = cell(length(sampling),1);
all_on_dot = zeros(length(sampling),length(fn));
all_dot_n = zeros(length(sampling),length(fn));
for i=1:length(fn)
    try
        load(fullfile(path,fn{i},'count_with_dot.mat'),'count_data');
        on_dot = [count_data.on_dot];
        clear count_data
        load(fullfile(path,fn{i},[fn{i} '_count']),'count_data');
        for j=1:length(count_data)
            count_data(j).on_dot = on_dot(j);
        end
    catch
        load(fullfile(path,fn{i},[fn{i} '_count']),'count_data');
    end
    load(fullfile(path,fn{i},fn{i}),'mosq');
    if strcmp(fn{i},'20171019-orl1010-0')
        [samplepos,on_dot,dot_n] = gather_pos(count_data,dotpos,polypos,sampling,interval,0);
    else
        [samplepos,on_dot,dot_n] = gather_pos(count_data,dotpos,polypos,sampling,interval,1);
    end
    dot_n = dot_n / (interval*mosq);
    
    for j=1:length(sampling)
        allpos{j} = vertcat(allpos{j},samplepos{j});
    end
    all_on_dot(:,i) = on_dot;
    all_dot_n(:,i) = dot_n;
end

end
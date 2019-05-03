

function [samplepos,on_dot,dot_n] = gather_pos(count_data,dotpos,polypos,sampling,interval,pos_counted)
% start: time point to start (heat bout is 540+90=630)
% interval: seconds to gather (bout is 90)
% assumed length of trial is 9180, assumed interbout interval is 720

if nargin < 6
    pos_counted = 1;
    if nargin < 5
        interval = 90;
        if nargin < 4
            sampling = 630:720:9180;
        end
    end
end

samplepos = cell(length(sampling),1);
on_dot = zeros(length(sampling),1);
dot_n = zeros(length(sampling),1);

for s=1:length(sampling)
    % collect listed positions
    pos = vertcat(count_data(sampling(s):sampling(s)+interval).pos);
    % find positions on dots and randomly assign positions by drawing from
    % dotpos
    dot_n(s) = sum([count_data(sampling(s):sampling(s)+interval).on_dot]);
    samplepos{s} = pos;
    counted_dot = sum(inpolygon(pos(:,1),pos(:,2),polypos(:,1),polypos(:,2)));
    if pos_counted && (counted_dot > 0)
        dot_n(s) = dot_n(s) + counted_dot;
        samplepos{s} = vertcat(pos,dotpos(randi(length(dotpos),1,dot_n(s)),:));
    end
    on_dot(s) = dot_n(s) / (sum([count_data(sampling(s):sampling(s)+interval).right_pelt]) + dot_n(s));
    
end

end
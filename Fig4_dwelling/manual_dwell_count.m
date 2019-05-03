function [time,landing_n,takeoff_n,a] = manual_dwell_count(a,path,delay)
% a is a segment of dwell

if nargin < 3
    delay = 0;
    if nargin < 2
        path = 'D:\Dropbox\labwork\behavior\heatseeking-pilot\';
    end
end

dir_name = a(1).dir_name;
cur_dir = [path dir_name];
%load(fullfile(cur_dir,[dir_name '_count']),'count_data');

start = [];
stop = [];
for i=(delay+1):length(a)
    [a(i).landing_n,~] = size(a(i).landing);
    [a(i).takeoff_n,~] = size(a(i).takeoff);

    %a(i).dmoz = diff([count_data(n-1:n).right])+diff([count_data(n-1:n).on_dot]);
    if i==delay+1
        a(i).dmoz = a(i).right;
    else
        a(i).dmoz = diff([a(i-1:i).right]);
    end
    dcount = a(i).landing_n - a(i).takeoff_n;
    if dcount < a(i).dmoz % underestimating landings
        a(i).landing_n = a(i).landing_n + (a(i).dmoz - dcount);
    elseif dcount > a(i).dmoz % underestimating takeoffs
        a(i).takeoff_n = a(i).takeoff_n + (dcount - a(i).dmoz);
    end
    
    start = [start repmat(a(i).fn,[1 a(i).landing_n])];
    stop = [stop repmat(a(i).fn,[1 a(i).takeoff_n])];

    if i==delay+1
        start = [start repmat(a(i).fn, [1 a(i).right + a(i).takeoff_n - a(i).landing_n])];
    elseif i==length(a)
        stop = [stop repmat(a(i).fn, [1 a(i).right])];
    end
end

if length(start) == length(stop)
    time = (sum(stop)-sum(start))/length(stop);
else
    disp('Values of landing and takeoff do not match up!');
    time = [];
end
landing_n = sum([a.landing_n]);
takeoff_n = sum([a.takeoff_n]);

plot([a(delay+1:end).fn],[a(delay+1:end).landing_n],'g-');
plot([a(delay+1:end).fn],[a(delay+1:end).takeoff_n],'m-');

end
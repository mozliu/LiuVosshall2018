function h = heatmap_wrapper(pos,path,name,ifdot)

for i=1:length(pos)
    if ifdot
        f = (pos{i}(:,1) == 298);
        newX = pos{i}(:,1);
        newY = pos{i}(:,2);
        pos{i} = horzcat(newX(~f),newY(~f));
    end
    h(i) = heatseeking_heatmap(pos{i},path,[name '-' num2str(i)],ifdot);
end

savefig(h,[path '/' name]);

end
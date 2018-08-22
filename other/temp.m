

temp = cellfun(@(x) squeeze(x{1}(end,1,:)), {kinData.modifiedLocations}, 'UniformOutput', false);
temp = cat(2, temp{:})';

figure;
for i = 1:size(temp,1)
    plot(temp(i,:)); hold on;
end
pimpFig
function plotSamples () 
col = (0:1/215:1)'*[1,1,1];
f = 0:3;
for i = 1:length(f)
    subplot(2,2,i);
    image(setFocus(f(i)));
    colormap(col);
    axis('equal');
    axis('off');
    title(sprintf('Image with focus set to %1.1f',f(i)));
end
print -depsc '../PICTURES/plotSamples.eps'

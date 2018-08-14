function colourFlip(h,cmap)

% function to flip the order of faceColors in an area graph, h being the
% handle for the area plot, cmap being the current colour map for the
% figure.


Lh = length(h);

Lc = length(cmap);

Id = [round((Lc/Lh)*(Lh:-1:1))];
    
for i = 1:Lh
     set(h(i),'FaceColor',cmap(Id(i),:));
end
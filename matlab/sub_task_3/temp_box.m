function box=temp_box(y,x,ylimit,xlimit,edge_length)
box=zeros(edge_length,edge_length,2,'double');
for yside=1:edge_length
    for xside=1:edge_length
        x_coor=x-(edge_length-1)/2+xside-1;
        y_coor=y-(edge_length-1)/2+yside-1;
        if x_coor<=0,xcoor=xlimit-x_coor;end
        if y_coor<=0,ycoor=ylimit-y_coor;end
        if x_coor>xlimit,xcoor=x_coor-xlimit;end
        if y_coor>xlimit,ycoor=y_coor-ylimit;end
        box(yside,xside,:)=[x_coor y_coor];
    end
end

end
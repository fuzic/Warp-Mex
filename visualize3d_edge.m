function []=visualize3d_edge(im)
    % im is a 4d BINARY edge-space object
    [az,el]=view();
    %figure;
    clf;
    grid on;
    hold on;
    xs=[];
    ys=[];
    zs=[];
    cs=[];
    cmap=colormap(rand(1+max(max(max(max(im)))),3));
    for x=1:size(im,1)
        for y=1:size(im,2)
            for z=1:size(im,3)
                if im(x,y,z,1)~=0
                    xs(end+1,:)=[x,x-1];
                    ys(end+1,:)=[y,y];
                    zs(end+1,:)=[z,z];
                    display(im(x,y,z,1))
                    cs(end+1,:)=cmap(ceil(im(x,y,z,1)+1),:);
                end
                if im(x,y,z,2)~=0
                    xs(end+1,:)=[x,x];
                    ys(end+1,:)=[y,y-1];
                    zs(end+1,:)=[z,z];
                    cs(end+1,:)=cmap(ceil(im(x,y,z,2)+1),:);
                end
                if im(x,y,z,3)~=0
                    xs(end+1,:)=[x,x];
                    ys(end+1,:)=[y,y];
                    zs(end+1,:)=[z,z-1];
                    cs(end+1,:)=cmap(ceil(im(x,y,z,3)+1),:);
                end
            end
        end
    end
    for i=1:size(xs,1)
	    plot3(xs(i,:),ys(i,:),zs(i,:),'Color',cs(i,:));
	end
    axis([0,size(im,1),0,size(im,2),0,size(im,3)]);
    view(az,el);
    hold off;

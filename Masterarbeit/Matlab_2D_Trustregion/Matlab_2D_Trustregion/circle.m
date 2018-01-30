function [  ] = circle( x,y,r )
%CIRCLE Summary of this function goes here
%   Detailed explanation goes here
rectangle('Position',[x-r,y-r,2*r,2*r],'Curvature',[1,1],'EdgeColor','g');

end


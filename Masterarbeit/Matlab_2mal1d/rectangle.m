function [  ] = rectangle( mu_x,mu_y,fn,kt )
%RECTANGLE Summary of this function goes here
%   Detailed explanation goes here
    hold on
    h1=line([-mu_x*fn/kt,mu_x*fn/kt],[mu_y*fn/kt,mu_y*fn/kt],'LineWidth',1,'LineStyle','--');
    set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');  

    hold on
    h2=line([-mu_x*fn/kt,-mu_x*fn/kt],[mu_y*fn/kt,-mu_y*fn/kt],'LineWidth',1,'LineStyle','--');
    set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');  

    hold on
    h3=line([mu_x*fn/kt,mu_x*fn/kt],[mu_y*fn/kt,-mu_y*fn/kt],'LineWidth',1,'LineStyle','--');
    set(get(get(h3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');  

    hold on
    h4=line([-mu_x*fn/kt,mu_x*fn/kt],[-mu_y*fn/kt,-mu_y*fn/kt],'LineWidth',1,'LineStyle','--');
    set(get(get(h4,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');  


end


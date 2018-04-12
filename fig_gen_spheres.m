function fig_gen_spheres(pos_mag_2,angle_view,color_arrow,shift_in_z)
if nargin<4
    shift_in_z=0;
end
pos_mag=pos_mag_2(1,:);
if size(pos_mag_2,1)>1
    pos_mag0=pos_mag_2(2,:);
else
    pos_mag0=[0 0 0];
end
pos_mag0(:,3)=pos_mag0(:,3)+shift_in_z;
pos_mag(:,3)=pos_mag(:,3)+shift_in_z;
if nargin<3
    color_arrow='k';
end

%axes('Position',[0.3 0.3 4 4]);
set(gcf,'Units','pixels');
pos = set(gcf,'Position',round([10 10 300 450]*2));

if nargin<2
    view([51 31]);
else
    view(angle_view);
end

inca=1;
list_angl=(0:inca:360)*pi/180;
hold on
plot3(0*cos(list_angl),1*cos(list_angl),1*sin(list_angl)+shift_in_z,'k-','color',[1 1 1]*0.5)
plot3(1*cos(list_angl),0*cos(list_angl),1*sin(list_angl)+shift_in_z,'k-','color',[1 1 1]*0.5)
plot3(1*cos(list_angl),1*sin(list_angl),0*cos(list_angl)+shift_in_z,'k-','color',[1 1 1]*0.5)
plot3([0 0],[0 0],[-1 1]+shift_in_z,'k-','color',[1 1 1]*0.5)
plot3([0 0],[-1 1],[0 0]+shift_in_z,'k-','color',[1 1 1]*0.5)
plot3([-1 1],[0 0],[0 0]+shift_in_z,'k-','color',[1 1 1]*0.5)
axis('equal')
%disp(['offset_first_null : ' num2str(offsset_first_null) ' Hz ' ])
if shift_in_z==0
for loarrows=size(pos_mag,1)
    
    
    %plot field vector
    
    % if disp_on
    % disp(['Offset : ' num2str(loop_offset) ' Hz w_eff=' num2str(nu_eff) ' Hz'])
    
    %   plot3([sin(tilt_angle)],[ 0],[ cos(tilt_angle)],'ko','MarkerSize',5)
    
    % end
    
    %%%%  plot3(stor_tr(:,1),stor_tr(:,2),stor_tr(:,3),'k:','linewidth',1.25)
    
    %  if disp_on,
    % plot3([pos_mag_p(1,1) ],[pos_mag_p(1,2) ],[ pos_mag(1,3)],'ko')
    % plot3([pos_mag_p(1,1) ],[pos_mag_p(1,2) ],[ pos_mag(1,3)],'k.','MarkerSize',12)
    plot3([0 pos_mag(loarrows,1) ],[0 pos_mag(loarrows,2) ],[0 pos_mag(loarrows,3)]+shift_in_z)
    plot3([0 pos_mag(loarrows,1) ],[0 pos_mag(loarrows,2) ],0*[0 pos_mag(loarrows,3)]+shift_in_z,'k:')
    plot3(pos_mag(loarrows,1) *[1 1],pos_mag(loarrows,2)*[1 1 ],1*[0 pos_mag(loarrows,3)]+shift_in_z,'k:')
  %  plot3([0 pos_mag(loarrows,1) ],[ pos_mag(loarrows,2) pos_mag(loarrows,2) ],0*[0 pos_mag(loarrows,3)],'g-')
  %  plot3([pos_mag(loarrows,1) pos_mag(loarrows,1) ],[0 pos_mag(loarrows,2) ],0*[0 pos_mag(loarrows,3)],'b-')
   % plot3([pos_mag(loarrows,1) pos_mag(loarrows,1) ],[pos_mag(loarrows,2) pos_mag(loarrows,2) ],[0 pos_mag(loarrows,3)],'r-')
    arrow(pos_mag0,[pos_mag(loarrows,1)  pos_mag(loarrows,2)  pos_mag(loarrows,3)],'linewidth',1.5,'color',color_arrow)
    %  end
    drawnow
    
    % plot3([start_pt_p(1,1) start_pt(1,1)],[start_pt_p(1,2) start_pt(1,2)],[start_pt_p(1,3) start_pt(1,3)],'g-')
    
end
end
end

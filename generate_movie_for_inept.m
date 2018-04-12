clear all
every_pt=500;
%every_pt=2500;% this is to speed up generation
for center_antiphase_arrow=[0 1]
for loop_j=0
    for mainlooop=0:6
        
        if  mainlooop==0
            text_title='Pulse Evolution PI(y) inv(J) Evolution (no chemical shift evolution)';
            file_name=['Spin_echo_for_INEPT_PIy_NO_chem_shift_evol.avi'];
            
        end
        if  mainlooop==1
            text_title='Pulse Evolution PI(x) inv(J) Evolution (no chemical shift evolution)';
            file_name=['Spin_echo_for_INEPT_PIx_NO_chem_shift_evol.avi'];
            
        end
        if  mainlooop==2
            text_title='Pulse Evolution PI(y) inv(J) Evolution ';
            file_name=['Spin_echo_for_INEPT_PIy_chem_shift_evol.avi'];
            
        end
        if  mainlooop==3
            text_title='Pulse Evolution PI(x) inv(J) Evolution ';
            file_name=['Spin_echo_for_INEPT_PIx_chem_shift_evol.avi'];
            
        end
        if  mainlooop==4
            text_title='Pulse Evolution PI(y) inv(J) Evolution (imperfect tau value)';
            file_name=['Spin_echo_for_INEPT_PIy_chem_shift_evol_imperfect_tau.avi'];
            
        end
        if  mainlooop==5
            text_title='Pulse Evolution PI(x) inv(J) Evolution (imperfect tau value)';
            file_name=['Spin_echo_for_INEPT_PIx_chem_shift_evol_imperfect_tau.avi'];
            
        end
        if every_pt==2500
            file_name=['delete_quick_' file_name];
        end
         if center_antiphase_arrow
            file_name=['Centered_antiphase_arrow_' file_name];
        end
        writerObj = VideoWriter([file_name]);
        
        % for mainlooop=-3
        clear stor_tr
        clear stor_tr_crude
        clear stor_t_crude
        clear stor_t
        main_ratio=90;
        
        figure(1)
        %  Z = peaks;
        %  surf(Z)
        %  axis tight manual
        %  ax = gca;
        %  ax.NextPlot = 'replaceChildren';
        
        
        %  writerObj = VideoWriter(['mov_' num2str(mainlooop) '.mp4']);
        writerObj.Quality = 75;%75 is default
        
        
        open(writerObj);
        
        axis tight
        set(gca,'nextplot','replacechildren');
        set(gcf,'Renderer','zbuffer');
        
        list_of_values= 0.00:1:main_ratio+1;
        loops = 360;
        loops = 1;
        count_sto=1;
        F(size(list_of_values,2)+1) = struct('cdata',[],'colormap',[]);
        counter_frames=1;
        
        pul_dur=10e-6;
        angle_pulse=90/180*pi;%90 deg
        ampli_hz=(angle_pulse/pul_dur)/(2*pi);
        disp(['pulse amplitude : ' num2str(ampli_hz) ' Hz'])
        loop_offset=0+000*1000;
        start_pt=[0 -1 0];
        evol_time=0.1;
        pulse_time=0.04;
        tmax=evol_time+pulse_time;
        center=10.123421;
        center=0;
        if mainlooop>1
            center=0.33;
        end
        J=0.2*2*pi;%pefect refocus...
        J=1;%pefect refocus...
        if mainlooop>=4
            J=1.1;
        end
        
        larmo=[center ]+[-J/2 J/2];
        
        starting_vector=[ 0 -1 0];
        pos_mag=(starting_vector'*ones(1,size(larmo,2)))';%starts on x...
        t_step=0.05;
        firstr1=0;
        tau=10;%(defaul no second pulse
        
        
        % if mod(round(t*20),6)==0
        R2=0;
        R1=0;
        before_180=1;
        for t = 0:0.000001:tmax
            % for t = 0:0.000001:0.05
            % X = sin(loopj*pi/10)*Z;
            %surf(X,Z)
            
            %  larmo=[10.123122 ];
            
            increment_tilt=pi/100000;
            inc=0;
            if (t>evol_time/2) && (t<(evol_time/2+pulse_time))
                % apply pi pulse
                if (mainlooop ==1) || (mainlooop ==3) || (mainlooop ==5)
                    %        pos_mag(:,1)=-pos_mag(:,1);before_180=0;
                    rfy=5/(pulse_time/0.02);
                rfx=0.0000000;
                else
                    %       pos_mag(:,2)=-pos_mag(:,2);before_180=0;
                    rfx=5/(pulse_time/0.02);
                rfy=0.0000000;
                end
                
                allow_evol=0;
                
            else
                rfx=0.0000000;
                
                rfy=0.0000000;
                allow_evol=1;
            end
           
            %  tilt_angle=atan((ampli_hz/loop_offset));
            
            %  if tilt_angle<0, tilt_angle=tilt_angle+pi;end
            
            for loop=1:size(larmo,2)
                v=[rfx rfy allow_evol*larmo(1,loop) ];
                di=cross(v,pos_mag(loop,:));
                pos_mag(loop,:)=pos_mag(loop,:)+di*increment_tilt;
                
            end
            %  di=di/norm(pos_mag);
            
            %  pos_mag(1,1)=pos_mag(1,1)*exp(-R2*increment_tilt);
            %  pos_mag(1,2)=pos_mag(1,2)*exp(-R2*increment_tilt);
            if R1~=0
                if firstr1==0
                    valref=pos_mag(1,3);diff=1-valref;
                    firstr1=1;
                    
                    
                end
                diff=diff*((exp(-R1*increment_tilt)));
                pos_mag(1,3)=1-diff;
            end
            rf=rfx+rfy;
            nu_eff=sqrt(loop_offset*loop_offset+ampli_hz*ampli_hz);
            
            %
            %      how_much_further=2;
            %     if size(factor,2)==1
            %         plot3(how_much_further*[0 sin(tilt_angle)],[0 0],how_much_further*[0 cos(tilt_angle)],'k--')
            %         text(how_much_further*[sin(tilt_angle)],[ 0],how_much_further*[ cos(tilt_angle)],'Beff')
            %
            %     end
            %     %%%%%%%%%%%%%%%%%
            %
            %
            
            
            
            %    figure('Units','inches', 'PaperPositionMode','auto', 'Position',[0 0 4 4]);
            
            
            
            % end
            if before_180
                if (t>evol_time/2)
                    before_180=0;
                    
                    J=-J;
                    larmo=[center ]+[-J/2 J/2];
                end
                
            end
            % end
            
            
            %work on figure
            if mod(round(t*1e6),every_pt)==0
                
                
                
                where_low=-3;
                disp(['time: ' num2str(t) '/' max(num2str(tmax)) ' ' num2str(rf)] )
                cur_f=figure(1);clf;hold on
                
                % black arrow
                how_much_higher=-2.2;
                for loop=1:size(larmo,2)
                    
                    fig_gen_spheres(pos_mag(loop,:),[51 31],'k',0)
                    fig_gen_spheres(pos_mag(loop,:),[51 31],'k',how_much_higher)
                end
                %midd arro
                [v1 v2]=view();
                cen_vect=pos_mag(1,:)*0.5+0.5*pos_mag(2,:);
                color_inphase='g';
                color_antiphase='r';
                mi_dist=0.2;
                if sum(sum((cen_vect).*(cen_vect)))>(mi_dist/2)*(mi_dist/2)
                    fig_gen_spheres(cen_vect,[v1 v2],color_inphase)
                else
                    plot3( [0 cen_vect(:,1)],[0 cen_vect(:,2)],[ 0 cen_vect(:,3)],color_inphase,'LineWidth',2)
                end
                %in first sphere
                if sum(sum(((pos_mag(1,:)-cen_vect).*(pos_mag(1,:)-cen_vect))))>(mi_dist/2)*(mi_dist/2)
                    %at end...
                    %    fig_gen_spheres([  pos_mag(1,:);cen_vect],[v1 v2],color_antiphase)
                    %    fig_gen_spheres([  pos_mag(2,:);cen_vect],[v1 v2],color_antiphase)
                    %centered
                    fig_gen_spheres([  pos_mag(1,:);cen_vect]-cen_vect*center_antiphase_arrow,[v1 v2],color_antiphase)
                    fig_gen_spheres([  pos_mag(2,:);cen_vect]-cen_vect*center_antiphase_arrow,[v1 v2],color_antiphase)
                else
                    plot3( pos_mag(:,1),pos_mag(:,2),pos_mag(:,3),color_antiphase,'LineWidth',2)
                end
                 %in first sphere
                if sum(sum(((pos_mag(1,:)-cen_vect).*(pos_mag(1,:)-cen_vect))))>(mi_dist/2)*(mi_dist/2)
                    %at end...
                    %    fig_gen_spheres([  pos_mag(1,:);cen_vect],[v1 v2],color_antiphase)
                    %    fig_gen_spheres([  pos_mag(2,:);cen_vect],[v1 v2],color_antiphase)
                    %centered
                    fig_gen_spheres([  pos_mag(1,:);cen_vect]-cen_vect*center_antiphase_arrow,[v1 v2],color_antiphase)
                    fig_gen_spheres([  pos_mag(2,:);cen_vect]-cen_vect*center_antiphase_arrow,[v1 v2],color_antiphase)
                else
                    plot3( pos_mag(:,1),pos_mag(:,2),pos_mag(:,3),color_antiphase,'LineWidth',2)
                end
                
                 %on second sphere
                 pos_mag_other=pos_mag;
                 pos_mag_other(:,3)=sqrt(power(pos_mag_other(:,1)-cen_vect(1,1),2)+power(pos_mag_other(:,2)-cen_vect(1,2),2)+power(pos_mag_other(:,3)-cen_vect(1,3),2));
                 pos_mag_other(:,2)=0*pos_mag_other(:,3);
                 pos_mag_other(:,1)=0*pos_mag_other(:,3);
                if sum(sum(((pos_mag(1,:)-cen_vect).*(pos_mag(1,:)-cen_vect))))>(mi_dist/2)*(mi_dist/2)
                    %at end...
                    %    fig_gen_spheres([  pos_mag(1,:);cen_vect],[v1 v2],color_antiphase)
                    %    fig_gen_spheres([  pos_mag(2,:);cen_vect],[v1 v2],color_antiphase)
                    %centered
                    fig_gen_spheres([  pos_mag_other(1,:);[0 0 0]]+[0 0 how_much_higher],[v1 v2],color_antiphase)
                    fig_gen_spheres([  -pos_mag_other(1,:);[0 0 0]]+[0 0 how_much_higher],[v1 v2],color_antiphase)
                else
                    plot3( pos_mag_other(:,1),pos_mag_other(:,2),pos_mag_other(:,3)+[ how_much_higher],color_antiphase,'LineWidth',2)
                end
                
                
                
                % store last point for full trajectory plot
                stor_tr_crude(count_sto,:)=pos_mag(1,:);
                stor_t_crude(count_sto,:)=t;
                count_sto=count_sto+1;
                % interpolation of
                [stor_tr, stor_t]=interp_to_smooth(stor_tr_crude,stor_t_crude);
                
                
                
                % plot3(stor_tr(:,1),stor_tr(:,2),stor_tr(:,3),'k-','linewidth',1.25)
                
%                 plot3( [0.5 0.5],[-1 1],[where_low where_low],'k-')
%                 plot3(-[0.5 0.5],[-1 1],[where_low where_low],'k-')
%                 plot3(-[0.0 0.0],[-1 1],[where_low where_low],'k:')
%                 plot3(+[1 1],[-1 1],[where_low where_low],'k:')
%                 plot3([-1 1],[-1 -1],[where_low where_low],'k-')
%                 plot3(-[1 1],[-1 1],[where_low where_low]+0.5,'k-')
%                 plot3(-[1 1],[-1 1],[where_low where_low]+1,'k:')
%                 plot3(-[1 1],[-1 1],[where_low where_low]+0,'k:')
%                 plot3(-[1 1],[-1 -1],[0 1 ]+where_low,'k-')
%                 
                rf_arrow=1.0;
                if rfx>0
                    plot3( [-rf_arrow rf_arrow],[0 0],[ 0 0],'c','LineWidth',2)
                    text( [ rf_arrow]*1.1,[0 ],[  0],'Inversion pulse')
                end
                if rfy>0
                    plot3( [0 0],[-rf_arrow rf_arrow],[ 0 0],'c','LineWidth',2)
                     text([0 ], [ rf_arrow]*1.1,[  0],'Inversion pulse')

                end
                %   list_t=2*(t/tmax)*(1:size(stor_tr,1))/size(stor_tr,1)-1;
                list_t=2*(t/tmax)*(stor_t)/max(abs(stor_t))-1;
                %      plot3( 0.5+0.5*stor_tr(:,1),list_t,stor_tr(:,3)*0+where_low,'g-','linewidth',1.25)
                
                
                %      plot3(-0.5+0.5*stor_tr(:,2),list_t,stor_tr(:,3)*0+where_low,'b-','linewidth',1.25)
                %     plot3(stor_tr(:,3)*0-1,list_t,0.5+0.5*stor_tr(:,3)+where_low,'r-','linewidth',1.25)
                
                % phase_a=cos(anglet/360*2*pi)+j*sin(anglet/360*2*pi);%as complex number
                %  phase_a=phase_a/(abs(phase_a));
                %  [dist_in_hz erro_in_deg]=shap_fn3d(loopj-1,mainlooop,main_ratio);
                %  [dist_in_hz erro_in_deg]=shap_fn3d(loopj-1,mainlooop,main_ratio);
                
                
                axis([ -1     1    -1     1    where_low-0.3     1])
                %plot3(store_traj(:,1),store_traj(:,2),store_traj(:,3),'k-','linewidth',1.5)
                %plot3(store_traj(:,1),store_traj(:,2),0*store_traj(:,3),'k-','linewidth',1)
                axis off
                
                set(gcf,'color','w');
                text(0,0,1.3,text_title,'HorizontalAlignment','center')
                
                drawnow;
                tmp_frame = getframe;
                if counter_frames==1
                    si=size(tmp_frame.cdata);
                else
                    tmp_frame.cdata=tmp_frame.cdata(1:si(1,1),1:si(1,2),1:si(1,3));
                end
                
                F(counter_frames)=tmp_frame;
                writeVideo(writerObj,F(counter_frames));
                counter_frames=counter_frames+1;
            end
            
        end
        
        close(writerObj);
        %movie(F,2)
        
        %     figure(111)
        %     clf
        %     plot(store_dis_in_hz(:,1),store_erro_in_deg(:,1),'b-');
        %     hold on
        %     plot(store_dis_in_hz(:,2),store_erro_in_deg(:,2),'r-');
        %     print('-depsc','-tiff','-r600',[ 'Phase_error_nearby_small_signals' num2str(main_ratio)  '.eps']);%here
        
    end
end
end
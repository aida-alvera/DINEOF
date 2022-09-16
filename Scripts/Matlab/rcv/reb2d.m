function [MOut,flag]=reb2d(coast,longini,dlong,latini,dlat,indicei,indicej,vectdat,m,n,pmax,notmattpred,writegher,writemumm,visu,filenameout1,filenameout2,title1,showinlog10,climdat,climin,climax,maskorir)

       % display('min i'), min(indiceii)
       % max(indiceii)
       % min(indicejj)
       % max(indicejj)
       % size(vectdat)
        sparseout=sparse(indicei,indicej,vectdat,m,n,pmax);

        matrixOut=full(sparseout);
        matrixOut(notmattpred)=NaN;
        MOut=matrixOut;
        %mindata=min(min(matrixOut))
        %maxdata=max(max(matrixOut))
        if writegher==1
            filenameout=[filenameout2,'.gher'];
            gwrite(filenameout,matrixOut);
        end
        if writemumm==1
            matrixOutm=matrixOut';
            filenameout=[filenameout2,'.bin'];
            fidD = fopen(filenameout,'w','ieee-le'); 
            count = fwrite(fidD,matrixOutm,'float32');
	        status = fclose(fidD);            
        end
        if visu==1
                matrixOutm=matrixOut';  % to cancel for quicker production when combined with wrtemumm=1
                matrixOut=imrotate((matrixOutm),90);
%                matrixOut=imrotate((matrixOut'),90);
                b=longini+dlong*[1:size(matrixOut,2)]; % moove out of subroutine as uniform for all the run reconstructions?
                a=(latini-dlat*size(matrixOut,1))+dlat*[1:size(matrixOut,1)]; %oove out as uniform for all the run reconstructions?
                figure
                pcolor(b,a,matrixOut);  
                title(title1,'FontWeight','bold','FontSize',16)
                xlabel('longitude ','FontWeight','bold','FontSize',16) %here special
                ylabel('latitude','FontWeight','bold','FontSize',16) %here special
                shading flat, colorbar; 
                if showinlog10==1
                    set(gca,'CLim',[-1,2]) 
                    set(colorbar,'YTick',[-1:1:2],'YTickLabel',[0.1;1;10;100],'FontWeight','bold','FontSize',16)
                else
                    if climdat==1
                        set(gca,'CLim',[climin,climax]) %here special
                    end
                    set(colorbar,'FontWeight','bold','FontSize',16)
                end
                set(gca,'FontWeight','bold','FontSize',16)
                hold on
                %contour(b,a,maskorir,[0.5,0.5],'-k')
                %min(a)
                %max(a)
                plot(coast(:,1),coast(:,2),'k','linewidth',0.5); %here special
                %colorbar;
%                formatimg='png';
%                saveas(gcf,filenameout1,formatimg)
                %print('-depsc','-loose','-r300',[fileout_error_data]);
                filenameout=[filenameout1,'.png'];
                print('-dpng','-r300',[filenameout]);
        end
        flag=1;
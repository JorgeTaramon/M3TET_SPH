function plot_tetra(FigNo,GCOORD,EL2NOD,el,flag,color,linestyle,linewidth)

figure(FigNo);
tetramesh(EL2NOD(1:4,el)',GCOORD','EdgeColor',color,'FaceColor','none',...
    'LineStyle',linestyle,'LineWidth',linewidth);
fs = 16;
sc = 0;
xlabel('X');ylabel('Y');zlabel('Z');
grid on
axis equal

elxyz = GCOORD(1:3,EL2NOD(:,el));

% scatter3(elxyz(1,:),elxyz(2,:),elxyz(3,:),50,'r')
if fs>0 && flag>0
    switch flag
        case 1
            for inod=1:4
                text(elxyz(1,inod),elxyz(2,inod),elxyz(3,inod),num2str(inod),...
                     'FontSize',fs,'FontWeight','bold','Color',color);
            end
        case 2
            text(sum(elxyz(1,1:4))/4,sum(elxyz(2,1:4))/4,sum(elxyz(3,1:4))/4,num2str(el),...
                 'FontSize',fs,'FontWeight','bold');
        case 3
            for inod=1:4
                text(elxyz(1,inod),elxyz(2,inod),elxyz(3,inod),num2str(EL2NOD(el,inod)),...
                     'FontSize',fs,'FontWeight','bold');
            end
        case 4
            for inod=1:10
                text(elxyz(1,inod),elxyz(2,inod),elxyz(3,inod),num2str(inod),... EL2NOD(el,inod)),...
                     'FontSize',fs,'FontWeight','bold','Color',color);
            end
        case 5
            nVnod = length(unique(EL2NOD(:,1:4)));
            for inod=5:10
                text(elxyz(1,inod),elxyz(2,inod),elxyz(3,inod),num2str(EL2NOD(el,inod)-nVnod),...
                     'FontSize',fs,'FontWeight','bold');
            end
        case 5
            for inod=4:10
                text(elxyz(1,inod),elxyz(2,inod),elxyz(3,inod),num2str(inod-4),...
                     'FontSize',fs,'FontWeight','bold');
            end
    end              
end

if sc>0
    scatter3(GCOORD(1,EL2NOD(el,1)),...
             GCOORD(2,EL2NOD(el,1)),...
             GCOORD(3,EL2NOD(el,1)),...
             sc,'k','filled')
end

disp(' ')

end
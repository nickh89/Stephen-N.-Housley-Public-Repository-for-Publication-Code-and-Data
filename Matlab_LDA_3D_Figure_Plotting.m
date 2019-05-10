%%%% 'Raw_data' should be in the format of m x n, where m is the samples (in
%%%% this case individual afferents and n is the feature space).
% Raw_Data=[];
% %%% 'labels' is a continuous variable discriminating group membership 
% labels =[];
% %%% 'canBase' is the nx3 reduced canonical basis derived from LDA analysis in rstudio code ('LDA Analysis') 
% canBase=[];
% %%% 'canData' is the true canonical variables i.e. reduced dimensionality
% %%% of the original Raw_Data
% canData=Raw_Data*canBase;


colors ={[0.501960	0.501960 0.501960], [0.243137255	0.36078	0.87843], [0.8588	0.015686	0.0823529], [0.5764705	0.23137	0.8274509]}; % matched colors to raw electrophysiological recordings
figure('pos',[10 10 1200 1500])

%%%% construct the 3D ellipses 
X=zeros(21,21,4); Y=X; Z=X;
for i1=1:4 
    [X(:,:,i1),Y(:,:,i1),Z(:,:,i1)]=errorellipse(canDataIa(labelsIa==i1,:));
end
for i1=1:4 
    patch(surf2patch(X(:,:,i1),Y(:,:,i1),Z(:,:,i1)), 'Facecolor', colors{i1});
end
set(gca,'color','none','FontName', 'Arial','FontWeight','bold','FontSize',12)

% xticks(0:50:250)
% yticks(0:20:100)
% zticks(-20:20:120)

grid on
% 
% % xz
% for i2=1:4
%     SX=size(X(:,:,i2));
% Y3 = 4*ones(SX); % A projection along Y-axis by making all Y-values to 3.
% hold on
% surf1a = surf(X(:,:,i2),Y3,Z(:,:,i2));
% set(surf1a,'FaceAlpha',0,'edgecolor', colors{i2},'linewidth', 1,'facecolor', colors{i2});
% end
% 
% %backwall zy
% for i3=1:4
%     SX=size(X(:,:,i3));
% X3 = -4*ones(SX); % A projection along Y-axis by making all Y-values to 3.
% hold on
% surf1a = surf(X3,Y(:,:,i3),Z(:,:,i3));
% set(surf1a,'FaceAlpha',0,'edgecolor', colors{i3},'linewidth', 1, 'facecolor', colors{i3});
% end 
% 
% %floor xy
% 
% for i4=1:4
%     SX=size(X(:,:,i4));
% Z3 = -3*ones(SX); % A projection along Y-axis by making all Y-values to 3.
% hold on
% surf1a = surf(X(:,:,i4),Y(:,:,i4),Z3);
% set(surf1a,'FaceAlpha',0,'edgecolor', colors{i4},'linewidth', 1, 'facecolor', colors{i4});
% end 

xlabel('')
ylabel('')
zlabel('')
xticks('')
yticks('')
zticks('')

grid on


%%%% set angle of view
az = 23;
el = 17;
view(az, el);


light
% legend('WT', 'WT + OX', 'Apc^{-/-}', 'Apc^{-/-} + OX', 'Location', 'SouthWestOutside')
% % set(legend,'color','none');
 saveas(gcf, 'IaCan');
 export_fig IaCan -m5 -transparent;
% % % 
% OptionZ.FrameRate=10;OptionZ.Duration=20;OptionZ.Periodic=true;
% CaptureFigVid([16,45;0,0;85,45], 'IaCan',OptionZ)


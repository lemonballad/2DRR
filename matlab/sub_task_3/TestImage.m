range=-500:500;
sx=350;
sy=15;
x0=200;
nuy=1/150;
y0=400;
x=exp(-(range-x0).^2/2/sx);
y=exp(-sin(2*pi*nuy*range).^1/sy).*exp(-abs(range)/y0);
[Y,X]=ndgrid(y,x);
F=X.*Y;

% [myFile myPath] = uiputfile('*.ppt','Create New PPT File'); % New ppt file. Modify if you want to append to existing. Remove if you want to use current active presentation (see exportPreso comment below)
% filename = fullfile(myPath,myFile);
pptserver = actxserver('PowerPoint.Application'); % create ppt actx server
% exportPreso = invoke(pptserver.Presentations,'Add'); % Open new presentation. To append to existing ppt, use 'open'. To append to currentlyactive ppt, use 
exportPreso = get(pptserver,'ActivePresentation');
slideCount = get(exportPreso.Slides,'Count'); % Get page count of preso
slideCount = int32(double(slideCount)+1); % find the 'next' preso page
newSlide = invoke(exportPreso.Slides,'Add',slideCount,11); % create new slide
slideWidth = get(exportPreso.PageSetup,'SlideWidth'); % get slide width
slideHeight = get(exportPreso.PageSetup,'SlideHeight'); % get slide height
set(newSlide.Shapes.Title.TextFrame.TextRange,'Text',' '); % set slide title to blank


surf(range,range,F,'edgecolor','none');xlim([150 250]);ylim([-400 200]);zlim([0.4 2]);grid off;
labfont=16;lwidax=2;set(gca,'fontsize',labfont,'linewidth',lwidax,'color','none','ZColor',[1 1 1],'YColor',[1 1 1],'XColor',[1 1 1]);axis square;
xlabel('X DATA');ylabel('Y DATA');title('Test Mountain Appearance','Color',[1 1 1]);
set(gca,'Color','none');
set(gcf,'Color','none');

oldPosition = get(gca,'position');
oldRenderer = get(gcf,'renderer');
oldGCFcolor = get(gcf,'color');
oldGCAcolor = get(gca,'color');
set(gca,'position',get(gca,'outerposition')-get(gca,'tightinset')*[-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]); % make plot image 'tight' and remove whitespace. Seems to behave differently in Office 2007
set(gcf,'Renderer','OpenGL'); % change renderer to openGL
% set(gcf,'color','none');
% set(gca,'color','none');
print(gcf,'-dmeta','-r1200'); % copy figure image to clipboard in meta format
set(gca,'position',oldPosition); % change figure properties back
set(gcf,'Renderer',oldRenderer); 
set(gcf,'color',oldGCFcolor);
set(gca,'color',oldGCAcolor);
imageObject = invoke(newSlide.Shapes,'paste'); % page image to the new slide
set(imageObject.PictureFormat,'TransparentBackground','msoTrue'); % Set picture to transparent
% set(imageObject.PictureFormat,'TransparencyColor','16777215'); % Set whites to be transparent
set(imageObject.Fill,'Transparency',1); % Fill is transparent
set(imageObject.Fill,'Visible','msoFalse'); % make fill invisible
set(imageObject,'width',0.8*slideWidth); % set size of pasted image
imageWidth = get(imageObject,'width');
imageHeight = get(imageObject,'height');
imageTop = 1.2*(slideHeight-imageHeight)/2;
imageLeft = (slideWidth-imageWidth)/2;
set(imageObject,'Top',imageTop,'Left',imageLeft); % set position of image in ppt
% invoke(exportPreso,'SaveAs',filename,1); % Assumes new file. Modify to use an existing file. Remove if pasting to currently active preso
% invoke(exportPreso,'Close'); % Remove if pasting to currently active preso
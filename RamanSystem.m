function varargout = RamanSystem(varargin)
clc

% RAMANSYSTEM MATLAB code for RamanSystem.fig
%      RAMANSYSTEM, by itself, creates a new RAMANSYSTEM or raises the existing
%      singleton*.
%
%      H = RAMANSYSTEM returns the handle to a new RAMANSYSTEM or the handle to
%      the existing singleton*.
%
%      RAMANSYSTEM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RAMANSYSTEM.M with the given input arguments.
%
%      RAMANSYSTEM('Property','Value',...) creates a new RAMANSYSTEM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before RamanSystem_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to RamanSystem_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help RamanSystem

% Last Modified by GUIDE v2.5 27-Apr-2018 11:11:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @RamanSystem_OpeningFcn, ...
                   'gui_OutputFcn',  @RamanSystem_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before RamanSystem is made visible.
function RamanSystem_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to RamanSystem (see VARARGIN)

% Choose default command line output for RamanSystem
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
 cla(handles.axes1,'reset') ;
  cla(handles.axes2,'reset') ;
   cla(handles.axes3,'reset') ;
% UIWAIT makes RamanSystem wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = RamanSystem_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global xs1
global specin
[filename,pathname]=...
uigetfile({'*.*';'*.txt'},'selectdata'); 
str=[pathname filename]; 
specin1=load(str);
% axis on
specin=specin1(:,2);
specin=specin';
xs1=specin1(:,1);
xs1=xs1';
str5=floor(xs1(1));
str4=floor(xs1(length(xs1)));
str6=ceil(abs((xs1(1)-xs1(2))));
set(handles.text4,'String',str4,'FontSize',10);   
set(handles.text5,'String',str5,'FontSize',10);   
set(handles.text6,'String',str6,'FontSize',10);   
axes(handles.axes1); %使用第一个axes

plot(xs1,specin);

% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%%%%%%去基线
 cla(handles.axes2,'reset') ;
global xs1
global specin
    spikespec=specin;
    s1=abs(sgolayfilt(spikespec,0,3));
    
    %通过残差分布计算阈值参数
    res=abs(spikespec-s1);
    stempp=spikespec;
     %标准差
    stdsg=0.7413*iqr(res);%标准差估计
    stempp=stempp.*(abs(res)<3*stdsg)+s1.*(abs(res)>3*stdsg);
    specin1=stempp;
    max_ind_yn=[false specin1(2:end-1)>specin1(3:end)&specin1(2:end-1)>specin1(1:end-2)&specin1(2:end-1)>2  false];

w1=1;
xsc1=length(xs1);
for i=1:xsc1 
    if max_ind_yn(i)==1
         aa(w1)=i;
         w1=w1+1;
    end
    
        
end
[mw, wsize]=size(aa);

for i=1:wsize-1
    if aa(i)-9<1
    aa(i)=10;
    end
    if aa(i)+9>xsc1
    aa(i)=xsc1-9;
    end
    w=aa(i);
    
    
    bb=specin1(w)/2;
    
    xnh1=[xs1(w-2) xs1(w-2) xs1(w+1) xs1(w+2)];
    specinnh1=[specin1(w-2) specin1(w-2) specin1(w+1) specin1(w+2)];
    pfcoef=polyfit(specinnh1,xnh1,1);
    xi=polyval(pfcoef,bb);

    if abs(xi-xs1(w))<=3 %w处是尖锋，再想一想
       
       pfcoef=polyfit(xs1(:,(w-9:w-5)),specin1(:,(w-9:w-5)),1);
       specin1(:,(w-4:w-1))=polyval(pfcoef,xs1(:,(w-4:w-1)));
       pfcoef=polyfit(xs1(:,(w+5:w+9)),specin1(:,(w+5:w+9)),1);
       specin1(:,(w+1:w+4))=polyval(pfcoef,xs1(:,(w+1:w+4)));
       xnh=[xs1(w-2) xs1(w-2) xs1(w+1) xs1(w+2)];
       specinnh=[specin1(w-2) specin1(w-2) specin1(w+1) specin1(w+2)];
       pfcoef=polyfit(xnh,specinnh,1);
       specin1(w)=polyval(pfcoef,xs1(w));
      
  
      
    end
        


end
%%%%%%
specinf=sgolayfilt(specin1,0,3);
coefsx1=cwt(specinf,20,'bior1.3');
coefsx1=abs(coefsx1);
FX1= gradient(specinf);
max_pos1=[false (FX1(2:end-1)>0&FX1(3:end)<0) false];
minxb1=[false coefsx1(2:end-1)<coefsx1(3:end)&coefsx1(2:end-1)<coefsx1(1:end-2) false];
ww1=1;

for i=11:length(xs1)-10

%  biaozhun1=[specinf((i-10):(i+10))];
 biaozhunz1=specinf(i);
%  biaozhun1m=mean(biaozhun1);
  
 if minxb1(i)==1&(max_pos1(i)==0&max_pos1(i+1)==0&max_pos1(i-1)==0)
 
         aaxb1(ww1)=i;
         ww1=ww1+1;
       
 end
    
        
end


xsc1=length(xs1);
aaxby1=[1 aaxb1 xsc1 ];
ff1=length(aaxby1);
for i=1:ff1-1
   xsfd1=[xs1(aaxby1(i)) xs1(aaxby1(i+1))];
   specoutfd1=[specinf(aaxby1(i)) specinf(aaxby1(i+1))];
   specoutfdk1=specinf(aaxby1(i):aaxby1(i+1));
   estblfd1=interp1(xsfd1,specoutfd1,xs1(aaxby1(i):aaxby1(i+1)),'linear');
   corrected_spectrumfd1=specoutfdk1-estblfd1;
   for j=1:1:length(corrected_spectrumfd1)
   if (corrected_spectrumfd1(j)<0)
       
%       corrected_spectrumfd1(j)=0;
       corrected_spectrumfd1(j)=abs(corrected_spectrumfd1(j));
   end
   end
   corrected_spectrumfdz1(aaxby1(i):aaxby1(i+1))=corrected_spectrumfd1;
  estblfdz1(aaxby1(i):aaxby1(i+1))=estblfd1;
end
axes(handles.axes2);
axis on
plot(xs1,corrected_spectrumfdz1);
hold on
plot(xs1,specin,'--g');
hold off

% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%%%%去尖峰
    cla(handles.axes2,'reset') ;
    global xs1
    global specin
    spikespec=specin;
    s1=abs(sgolayfilt(spikespec,0,3));
    
    %通过残差分布计算阈值参数
    res=abs(spikespec-s1);
    stempp=spikespec;
     %标准差
    stdsg=0.7413*iqr(res);%标准差估计
    stempp=stempp.*(abs(res)<3*stdsg)+s1.*(abs(res)>3*stdsg);
    specin1=stempp;
    max_ind_yn=[false specin1(2:end-1)>specin1(3:end)&specin1(2:end-1)>specin1(1:end-2)&specin1(2:end-1)>2  false];

w1=1;
xsc1=length(xs1);
for i=1:xsc1 
    if max_ind_yn(i)==1
         aa(w1)=i;
         w1=w1+1;
    end
    
        
end
[mw, wsize]=size(aa);

for i=1:wsize-1
    if aa(i)-9<1
    aa(i)=10;
    end
    if aa(i)+9>xsc1
    aa(i)=xsc1-9;
    end
    w=aa(i);
    
    
    bb=specin1(w)/2;
    
    xnh1=[xs1(w-2) xs1(w-2) xs1(w+1) xs1(w+2)];
    specinnh1=[specin1(w-2) specin1(w-2) specin1(w+1) specin1(w+2)];
    pfcoef=polyfit(specinnh1,xnh1,1);
    xi=polyval(pfcoef,bb);

    if abs(xi-xs1(w))<=3 %w处是尖锋，再想一想
       
       pfcoef=polyfit(xs1(:,(w-9:w-5)),specin1(:,(w-9:w-5)),1);
       specin1(:,(w-4:w-1))=polyval(pfcoef,xs1(:,(w-4:w-1)));
       pfcoef=polyfit(xs1(:,(w+5:w+9)),specin1(:,(w+5:w+9)),1);
       specin1(:,(w+1:w+4))=polyval(pfcoef,xs1(:,(w+1:w+4)));
       xnh=[xs1(w-2) xs1(w-2) xs1(w+1) xs1(w+2)];
       specinnh=[specin1(w-2) specin1(w-2) specin1(w+1) specin1(w+2)];
       pfcoef=polyfit(xnh,specinnh,1);
       specin1(w)=polyval(pfcoef,xs1(w));
      
  
      
    end
        


end


axes(handles.axes2);
axis on
plot(xs1,specin1);

% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc
clear all


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 cla(handles.axes3,'reset') ;
global xs1
global specin
   spikespec=specin;
    s1=abs(sgolayfilt(spikespec,0,3));
    
    %通过残差分布计算阈值参数
    res=abs(spikespec-s1);
    stempp=spikespec;
     %标准差
    stdsg=0.7413*iqr(res);%标准差估计
    stempp=stempp.*(abs(res)<3*stdsg)+s1.*(abs(res)>3*stdsg);
    specin1=stempp;
    max_ind_yn=[false specin1(2:end-1)>specin1(3:end)&specin1(2:end-1)>specin1(1:end-2)&specin1(2:end-1)>2  false];

w1=1;
xsc1=length(xs1);
for i=1:xsc1 
    if max_ind_yn(i)==1
         aa(w1)=i;
         w1=w1+1;
    end
    
        
end
[mw, wsize]=size(aa);

for i=1:wsize-1
    if aa(i)-9<1
    aa(i)=10;
    end
    if aa(i)+9>xsc1
    aa(i)=xsc1-9;
    end
    w=aa(i);
    
    
    bb=specin1(w)/2;
    
    xnh1=[xs1(w-2) xs1(w-2) xs1(w+1) xs1(w+2)];
    specinnh1=[specin1(w-2) specin1(w-2) specin1(w+1) specin1(w+2)];
    pfcoef=polyfit(specinnh1,xnh1,1);
    xi=polyval(pfcoef,bb);

    if abs(xi-xs1(w))<=3 %w处是尖锋，再想一想
       
       pfcoef=polyfit(xs1(:,(w-9:w-5)),specin1(:,(w-9:w-5)),1);
       specin1(:,(w-4:w-1))=polyval(pfcoef,xs1(:,(w-4:w-1)));
       pfcoef=polyfit(xs1(:,(w+5:w+9)),specin1(:,(w+5:w+9)),1);
       specin1(:,(w+1:w+4))=polyval(pfcoef,xs1(:,(w+1:w+4)));
       xnh=[xs1(w-2) xs1(w-2) xs1(w+1) xs1(w+2)];
       specinnh=[specin1(w-2) specin1(w-2) specin1(w+1) specin1(w+2)];
       pfcoef=polyfit(xnh,specinnh,1);
       specin1(w)=polyval(pfcoef,xs1(w));
      
  
      
    end
        


end
%%%%%%
specinf=sgolayfilt(specin1,0,3);
coefsx1=cwt(specinf,20,'bior1.3');
coefsx1=abs(coefsx1);
FX1= gradient(specinf);
max_pos1=[false (FX1(2:end-1)>0&FX1(3:end)<0) false];
minxb1=[false coefsx1(2:end-1)<coefsx1(3:end)&coefsx1(2:end-1)<coefsx1(1:end-2) false];
ww1=1;

for i=11:length(xs1)-10

%  biaozhun1=[specinf((i-10):(i+10))];
 biaozhunz1=specinf(i);
%  biaozhun1m=mean(biaozhun1);
  
 if minxb1(i)==1&(max_pos1(i)==0&max_pos1(i+1)==0&max_pos1(i-1)==0)
 
         aaxb1(ww1)=i;
         ww1=ww1+1;
       
 end
    
        
end


xsc1=length(xs1);
aaxby1=[1 aaxb1 xsc1 ];
ff1=length(aaxby1);
for i=1:ff1-1
   xsfd1=[xs1(aaxby1(i)) xs1(aaxby1(i+1))];
   specoutfd1=[specinf(aaxby1(i)) specinf(aaxby1(i+1))];
   specoutfdk1=specinf(aaxby1(i):aaxby1(i+1));
   estblfd1=interp1(xsfd1,specoutfd1,xs1(aaxby1(i):aaxby1(i+1)),'linear');
   corrected_spectrumfd1=specoutfdk1-estblfd1;
   for j=1:1:length(corrected_spectrumfd1)
   if (corrected_spectrumfd1(j)<0)
       
      corrected_spectrumfd1(j)=abs(corrected_spectrumfd1(j)); 
   end
   end
   corrected_spectrumfdz1(aaxby1(i):aaxby1(i+1))=corrected_spectrumfd1;
  estblfdz1(aaxby1(i):aaxby1(i+1))=estblfd1;
end


specin1=corrected_spectrumfdz1;
specinf=sgolayfilt(specin1,0,3);
coefsx1=cwt(specinf,5,'db1');
coefsx1=abs(coefsx1);

FX1 = gradient(specinf);
tidu=[false (FX1(2:end-1)>0&FX1(3:end)<0)|(FX1(2:end-1)<0&FX1(3:end)>0) false];
minxb1=[false coefsx1(2:end-1)<coefsx1(3:end)&coefsx1(2:end-1)<coefsx1(1:end-2) false];
maxb1=[false specinf(2:end-1)>specinf(3:end)&specinf(2:end-1)>specinf(1:end-2) false];
ww1=1;
w1=1;
meanzong=mean(specinf);
for i=6:length(xs1)-6
    biaozhunz1=specinf(i);
    if (tidu(i)==1|tidu(i+1)==1|tidu(i-1)==1)&minxb1(i)==1&maxb1(i)==1&biaozhunz1>2*meanzong
    aaxb(w1)=i;
     w1=w1+1;   
         bb=specinf(i)/2; 
         pfcoef=polyfit(specinf(:,(i-5:i)),xs1(:,(i-5:i)),1);
         xi=polyval(pfcoef,bb);
         for j=1:1:i
            minn(j)= abs(xs1(j)-xi);
         end
         minn1=min(min(minn));
         indx=find(minn==minn1);
         indx1=abs(indx-i);
        
 
         biaozhun1=abs(2*specinf(i)-specinf(max(1,i-2*indx1))-specinf(min(length(xs1),i+2*indx1)));

         if (tidu(i)==1|tidu(i+1)==1|tidu(i-1)==1)&minxb1(i)&maxb1(i)==1&biaozhunz1<=2*biaozhun1
                 
            aaxb1(ww1)=i;
             ww1=ww1+1;
         end
         
     
       
    end
     
end


axes(handles.axes3);
% plot(xs1,specin,'--g');
axis on
hold on
plot(xs1(aaxb),specinf(aaxb),'r*');
plot(xs1,specinf);
hold off
pufeng1=xs1(aaxb);
save unknown.txt pufeng1 -ascii



% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%%%%%%%%%%%%%pca
cla(handles.axes1,'reset') ;
  cla(handles.axes2,'reset') ;
Mysample=[2.5 2.4;0.5 0.7;2.2 2.9;1.9 2.2;3.1 3.0;2.3 2.7;2 1.6;1 1.1;1.5 1.6;1.1 0.9];
Mysamplem=mean(Mysample);
[m n]=size(Mysample);
linemean=mean(Mysample);
%对所有样本进行处理
for i=1:n
    DataAdjust(:,i)=Mysample(:,i)-linemean(i);
end
%求特征协方差矩阵
covsample=cov( DataAdjust);
%求协方差矩阵的特征值和特征矩阵
[Vectors,Values] = eig(covsample);
for i=1:length(Values)
   Values1= max(Values');
end
%将特征值从大到小排列，选择其中K个，将其对应的特征向量组成特征向量矩阵
[Valuesp,iValuesp]=sort(Values1,'descend');
% [c,I]=max(Values1);
Vectorsp=Vectors(:,iValuesp);
j=1;
for i=1:length(Valuesp)
   accumulation(i)=sum(Valuesp(i))/sum(Valuesp);
   if accumulation(i)>=0.85
       rownumber(j)=i;
       j=j+1;
   end
end
 rownumber=[2 1];
for i=1:length(rownumber)
   Vectorspx(:,rownumber(i))=Vectorsp(:,rownumber(i));
 
end

%将样本点投影到选取的特征向量上
Finaldata=DataAdjust*Vectorspx;
%计算结束
axes(handles.axes1);
plot(Finaldata(:,1),Finaldata(:,2),'+');
hold on
k1=Vectors(1,1)./Vectors(2,1);
plot([-2,2],[-2*k1,2*k1]);
k2=Vectors(1,2)./Vectors(2,2);
plot([-2,2],[-2*k2,2*k2]);

xx = Finaldata(:,1); yy = Finaldata(:,2);
[x,y] = pol2cart(xx,yy);
k = convhull(x,y);
axes(handles.axes2);
plot(x(k),y(k),'r-',x,y,'b+')

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%%%%%pls
cla(handles.axes1,'reset') ;
  cla(handles.axes2,'reset') ;
load pz.txt   
mu=mean(pz);sig=std(pz);
rr=corrcoef(pz);   %求相关系数矩阵 
data=zscore(pz); 
n=3;m=3;   
x0=pz(:,1:n);y0=pz(:,n+1:end);  
e0=data(:,1:n);f0=data(:,n+1:end);  
num=size(e0,1);
chg=eye(n); 
for i=1:n 
    matrix=e0'*f0*f0'*e0; 
    [vec,val]=eig(matrix); %求特征值和特征向量 
    val=diag(val);  
    [val,ind]=sort(val,'descend'); 
    w(:,i)=vec(:,ind(1));    %提出最大特征值对应的特征向量 
    w_star(:,i)=chg*w(:,i);  
    t(:,i)=e0*w(:,i);       
    alpha=e0'*t(:,i)/(t(:,i)'*t(:,i)); 
    chg=chg*(eye(n)-w(:,i)*alpha');    
    e=e0-t(:,i)*alpha';     
    e0=e;    
    beta=t\f0;  %求回归方程的系数
    cancha=f0-t*beta;    
    ss(i)=sum(sum(cancha.^2)); 
    for j=1:num 
        t1=t(:,1:i);f1=f0; 
        she_t=t1(j,:);she_f=f1(j,:);  
        t1(j,:)=[];f1(j,:)=[];        
        beta1=[t1,ones(num-1,1)]\f1;  
        cancha=she_f-she_t*beta1(1:end-1,:)-beta1(end,:);  
        press_i(j)=sum(cancha.^2); 
    end 
    press(i)=sum(press_i); 
    Q_h2(1)=1; 
    if i>1, Q_h2(i)=1-press(i)/ss(i-1); end 
    if Q_h2(i)<0.0975 
        fprintf('提出的成分个数 r=%d',i); break 
    end 
end 
beta_z=t\f0;  
xishu=w_star*beta_z;   
mu_x=mu(1:n);mu_y=mu(n+1:end);  
sig_x=sig(1:n);sig_y=sig(n+1:end); 
ch0=mu_y-(mu_x./sig_x*xishu).*sig_y; 
for i=1:m 
    xish(:,i)=xishu(:,i)./sig_x'*sig_y(i);  
end 
sol=[ch0;xish] 
save mydata x0 y0 num xishu ch0 xish 
 load mydata 
ch0=repmat(ch0,num,1);    
yhat=ch0+x0*xish; 
y1max=max(yhat);
y2max=max(y0);
ymax=max([y1max;y2max])  
cancha=yhat-y0; 
% axes(handles.axes1);
% subplot(2,2,1) ;
% plot(0:ymax(1),0:ymax(1),yhat(:,1),y0(:,1),'*') 
% subplot(2,2,2) 
% plot(0:ymax(2),0:ymax(2),yhat(:,2),y0(:,2),'O') 
% subplot(2,2,3) 
% plot(0:ymax(3),0:ymax(3),yhat(:,3),y0(:,3),'H') 
axes(handles.axes2);
bar(xishu') 



% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%%%%%%%%%%lda
cla(handles.axes1,'reset') ;
  cla(handles.axes2,'reset') ;
cls1_data=[2.93 6.634; 2.53 7.79; 3.57 5.65;3.16 5.47];%第一个类的训练集
cls2_data=[2.58 4.44; 2.16 6.22; 3.27 3.52];%第二个类的训练集
%求期望
E_cls1=mean(cls1_data);%第一类数据的期望矩阵
E_cls2=mean(cls2_data);%第二类数据的期望矩阵
E_all=mean([cls1_data;cls2_data]);%所有训练集的期望矩阵
%%%%%%%%%%%%%%%%%%%%分类前画图%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(handles.axes1);
for i=1:4
     plot(cls1_data(i,1),cls1_data(i,2),'.r');
     hold on;
end;
plot(E_cls1(1),E_cls1(2),'^r');
hold on;
for i=1:3
     plot(cls2_data(i,1),cls2_data(i,2),'*b');
     hold on;
end;
plot(E_cls2(1),E_cls2(2),'^b');
plot(E_all(1),E_all(2),'vc');
hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%计算类间离散度矩阵：
x1=E_all-E_cls1;
x2=E_all-E_cls2;
Sb=4*x1'*x1/7+3*x2'*x2/7;%%%%%%%%%为什么不用式（3）？？？
%计算类内离散度矩阵
y1=0;
for i=1:4
    y1=y1+(cls1_data(i,:)-E_cls1)'*(cls1_data(i,:)-E_cls1);
end;
y2=0;
for i=1:3
    y2=y2+(cls2_data(i,:)-E_cls2)'*(cls2_data(i,:)-E_cls2);
end;
Sw=4*y1/7+3*y2/7;%%%%%%%%%为什么不用式（3）？？？
%求最大特征值和特征向量
[V,L]=eig(inv(Sw)*Sb);

% for i=1:length(L)
%    Values1= max(L');
% end
% [Valuesp,iValuesp]=sort(Values1,'descend');
% % [c,I]=max(Values1);
% Vectorsp=Vectors(:,iValuesp);



[as,bs]=min(min(L));
newspaces=V(:,bs);%最大特征值所对应的特征向量
ks=newspaces(2)/newspaces(1);
bs=E_all(2)-ks*E_all(1);
% plot([E_all(1) E_all(2)],[E_all(1)*ks+bs E_all(1)*ks+bs],'-b');%两点画出一条直线，画出较小特征值对应的特征向量
plot([4 2],[4*ks+bs 2*ks+bs],'-b');%4,2自己给出，为了画图方便
[a,b]=max(max(L));
newspace=V(:,b);%最大特征值所对应的特征向量
new_cls1_data=cls1_data*newspace;%训练后的数据集
new_cls2_data=cls2_data*newspace;%训练后的数据集
%%%%%%%%%%%%%%%%%%画图代码%%%%%%%%%%%%%%%%%

hold off
k=newspace(2)/newspace(1);
axes(handles.axes2);
plot([0,6],[0*k,6*k],'-c');%画出最大特征值对应的特征向量，即样本所组成的线性空间所投影的子空间
hold on
%0,6自己设定
% plot([E_all(1),E_all(2)],[0,6*ks],'-c');
% axis([0 7 -2 9]);
axis([ 2 6 0 4]);%自己设定防止坐标系，画图好看
%画出样本投影到子空间点
for i=1:4
    temp=cls1_data(i,:);
    newx=(temp(1)+k*temp(2))/(k*k+1);
    newy=k*newx;
    plot(newx,newy,'*r');
end;

for i=1:3
    temp=cls2_data(i,:);
    newx=(temp(1)+k*temp(2))/(k*k+1);
    newy=k*newx;
    plot(newx,newy,'ob');
end;
%预测
prediction=[4.81 3.46];
result=prediction*newspace;
temp=new_cls1_data-[result;result;result;result];
temp2=new_cls2_data-[result;result;result];
if(min(abs(temp))>min(abs(temp2)))
    output='该样本属于不合格产品'
else
    output='该样本属于合格产品'
end;



% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
global p5
val = get(hObject,'Value');
switch val
case 1
    p5=1;
case 2
    p5=2;
case 3
    p5=3;
case 4
    p5=4;
case 5
    p5=5;
end



 


% --- Executes during object creation, after setting all properties.
function uipanel5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes when selected object is changed in uipanel5.
function uipanel5_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel5 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
global p1

switch get(hObject,'Tag')
    case 'radiobutton3'
            p1='rigrsure';
       case 'radiobutton4'
            p1='heursure';
          case 'radiobutton1'
            p1='sqtwolog';
            case 'radiobutton2'
                 p1='minimaxi';
                 
end


% --- Executes when selected object is changed in uipanel3.
function uipanel3_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel3 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
global p2

switch get(hObject,'Tag')
    case 'radiobutton9'
            p2='s';
       case 'radiobutton8'
            p2='h';
end
         


% --- Executes when selected object is changed in uipanel4.
function uipanel4_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel4 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
global p3

switch get(hObject,'Tag')
    
    case 'radiobutton6'
            p3='one';
       case 'radiobutton7'
            p3='sln';
          case 'radiobutton5'
            p3='mln';
            
end
      


% --- Executes when selected object is changed in uipanel2.
function uipanel2_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel2 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
global p4

switch get(hObject,'Tag')
    case 'radiobutton14'
       prompt={'输入小波基:'};
       p=inputdlg(prompt,'input');  
       p4=strcat('db',p{1})
%        p4=str2num(p4);
            
       case 'radiobutton15'
          prompt={'输入小波基:'};
          p=inputdlg(prompt,'input');  
           p4=strcat('coif',p{1})
%            p4=str2num(p4);
           
          case 'radiobutton11'
               prompt={'输入小波基:'};
               p=inputdlg(prompt,'input');   
               p4=strcat('sym',p{1})
%                p4=str2num(p4);
            
            case 'radiobutton13'
                p4='dmey';
                
               case 'radiobutton12'
                  prompt={'输入小波基:'};
                  p=inputdlg(prompt,'input');  
                   p4=strcat('bior',p{1})
%                    p4=str2num(p4);
                 
                 case 'radiobutton10'
                   prompt={'输入小波基:'};
                   p=inputdlg(prompt,'input');  
                   p4=strcat('rbio',p{1})
%                    p4=str2num(p4);

                  
end
         


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
unk=load ('unknown.txt');
unk=unk';
% syms jusifuyixi
% fp=fopen('xxxxx.txt');
yizhi=importdata('xxxxx.txt');
% yizhi=yizhi';
yizhidata= yizhi.data;
str=yizhi.textdata;
nn1=length(unk);
[M,N]=size(yizhidata);

nzong=zeros(1,N);
for k=1:N
for i=1:nn1
     
    for j=1:M
       
       if abs(unk(i)-yizhidata(j,k))<=1    
           nzong(1,k)=nzong(1,k)+1;
       end
      
  
    end
  

end

end
[C,I] = max(nzong);
str=str(I);
if C>=5
 set(handles.text3,'String',str,'FontSize',10);   
else
    str1='数据库中不存在该种物质';
  set(handles.text3,'String',str1,'FontSize',10); 
end



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cla(handles.axes2,'reset') ;
global xs1
    global specin
  [c,s]=wavedec(specin,4,'sym8');

ca1=appcoef(c,s,'sym8',1);
cd1=detcoef(c,s,1);
sigma1=wnoisest(c,s,1);
n1=numel(cd1);
thr1=sigma1*sqrt(2*log(n1));
% thr1= wbmpen(c,s,sigma1,2)
% denoise1=wdencmp('gbl',c,s,'sym8',1,thr1,'s',1);
for i=1:1:n1
if abs(cd1(i))<=thr1
    cd1(i)=0;
else
%     cd1(i)=cd1(i);
 cd1(i)=(cd1(i)./abs(cd1(i)))*(abs(cd1(i))-thr1);
    
end
end

ca2=appcoef(c,s,'sym8',2);
cd2=detcoef(c,s,2);
sigma2=wnoisest(c,s,2);
% thr2= wbmpen(c,s,sigma2,2)
n2=numel(cd2);
thr2=sigma2*sqrt(2*log(n2));
% denoise2=wdencmp('gbl',c,s,'sym8',2,thr2,'s',1);
for i=1:1:n2
if abs(cd2(i))<=thr2
    cd2(i)=0;
else
%     cd2(i)=cd2(i);
 cd2(i)=(cd2(i)./abs(cd2(i)))*(abs(cd2(i))-thr2);
end
end

ca3=appcoef(c,s,'sym8',3);
cd3=detcoef(c,s,3);
sigma3=wnoisest(c,s,3);
% thr3= wbmpen(c,s,sigma3,2)
n3=numel(cd3);
thr3=sigma3*sqrt(2*log(n3));
% denoise3=wdencmp('gbl',c,s,'sym8',3,thr3,'s',1);
for i=1:1:n3
if abs(cd3(i))<=thr3
    cd3(i)=0;
else
%     cd3(i)=cd3(i);
 cd3(i)=(cd3(i)./abs(cd3(i)))*(abs(cd3(i))-thr3);
    
end
end


ca4=appcoef(c,s,'sym8',4);
cd4=detcoef(c,s,4);
sigma4=wnoisest(c,s,4);
% thr4= wbmpen(c,s,sigma4,2)
n4=numel(cd4);
thr4=sigma4*sqrt(2*log(n4));
% denoise4=wdencmp('gbl',c,s,'sym8',4,thr4,'s',1);
for i=1:1:n4
if abs(cd4(i))<=thr4
    cd4(i)=0;
else
%     cd4(i)=cd4(i);
 cd4(i)=(cd4(i)./abs(cd4(i)))*(abs(cd4(i))-thr4);
    
end
end
cx=[ca4,cd4,cd3,cd2,cd1];

specinf= waverec(cx,s,'sym8');
axes(handles.axes2);
axis on
plot(xs1,specinf);


% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
xlswrite('E:\matlab pro\gui\data0.xls',shuju);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text3.
function text3_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to text3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on key press with focus on pushbutton11 and none of its controls.
function pushbutton11_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

function varargout = AntColonyAlgorithmWithEdgeDetectionGui(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @AntColonyAlgorithmWithEdgeDetectionGui_OpeningFcn, ...
                   'gui_OutputFcn',  @AntColonyAlgorithmWithEdgeDetectionGui_OutputFcn, ...
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


function AntColonyAlgorithmWithEdgeDetectionGui_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);

function varargout = AntColonyAlgorithmWithEdgeDetectionGui_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;



function edit1_Callback(hObject, eventdata, handles)

function edit1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit2_Callback(hObject, eventdata, handles)

function edit2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pushbutton1_Callback(hObject, eventdata, handles)
%% Goruntuyu okutuyoruz ve double'a ceviriyoruz. Ardindan goruntunun satir
% ve sutun degerlerini belirliyoruz. 
global goruntu dosyaAdi
dosyaAdi = imgetfile;
goruntu = double(imread(dosyaAdi))./255;
h = handles.axes1;
axes(h);
imshow(goruntu);

function pushbutton2_Callback(hObject, eventdata, handles)
   
msgbox('Secilen image icin algoritma calistiriliyor... Lutfen bekleyin... ','Dikkat !','warn');
  
global goruntu dosyaAdi alpha beta phi rho
  
[goruntuSatir, goruntuSutun] = size(goruntu);
ToplamIterasyonSayisi=str2num(get(handles.edit1, 'String'));

tic;
%% gorunurluk islemine baslatiliyor. 
  
     goruntuKopyasi = zeros(size(goruntu));
     goruntuKopyasi_normal = 0;    
     
     for satir =1:goruntuSatir
        for sutun=1:goruntuSutun           
            %% temsili bir karinca takim olusturuluyor. bunun icin gerekli
            %% olan satir ve sutunlar belirleniyor. Filtrede denilebilir.
            taslak1 = [satir-2 sutun-1; satir-2 sutun+1; satir-1 sutun-2;
                       satir-1 sutun-1; satir-1 sutun; satir-1 sutun+1; 
                       satir-1 sutun+2; satir sutun-1];
            taslak2 = [satir+2 sutun+1; satir+2 sutun-1; satir+1 sutun+2; 
                       satir+1 sutun+1; satir+1 sutun; satir+1 sutun-1; 
                       satir+1 sutun-2; satir sutun+1];
            taslak0 = find(taslak1(:,1)>=1 & taslak1(:,1)<=goruntuSatir & taslak1(:,2)>=1 & taslak1(:,2)<=goruntuSutun & taslak2(:,1)>=1 & taslak2(:,1)<=goruntuSatir & taslak2(:,2)>=1 & taslak2(:,2)<=goruntuSutun);  
            % find fonksiyonu ile sifirdan farkli elemanlari tespit ediyoruz.
            taslak11 = taslak1(taslak0, :);
            taslak22 = taslak2(taslak0, :);
            taslak00 = zeros(size(taslak11,1));  % sifir matrisi olusturuyoruz.          
            
            for i = 1:size(taslak11,1)
                taslak00(i) = abs(goruntu(taslak11(i,1), taslak11(i,2)) - goruntu(taslak22(i,1), taslak22(i,2)));    % abs ile mutlak deger aliyoruz.
            end
            
            if size(taslak11,1) == 0
                goruntuKopyasi(satir, sutun) = 0;
                goruntuKopyasi_normal = goruntuKopyasi_normal + goruntuKopyasi(satir, sutun);
            else 
                sigma = 10;   
                            taslak00 = sin(pi .* taslak00./2./sigma);  
                goruntuKopyasi(satir, sutun) = sum(sum(taslak00.^2));
                goruntuKopyasi_normal = goruntuKopyasi_normal + goruntuKopyasi(satir, sutun);
            end
        end
     end 
     
%% normalizasyon islemi gerceklestiriliyor.
     goruntuKopyasi = goruntuKopyasi./goruntuKopyasi_normal; 
     goruntuKopyasi = goruntuKopyasi.*100;
     
%% feromon islemi icin feromon degeri belirleniyor.
     feromon = 0.0001 .* ones(size(goruntu));     % bir matrisi olusturuyoruz.    
%0.0001
%% olasilik formulu icin gerekli olan degerler belirtiliyor     
%      alpha = 1;   % alfa degeri belirleniyor.
%      beta = 1;    % beta degeri belirleniyor. 
%      rho = 0.1;   % rho degeri belirleniyor.
%      phi = 0.05;  % phi degeri belirleniyor.

%% goruntu uzerindeki satir sutun degerleri ile toplam karinca sayisi
% belirliyoruz. Her bir karinca icin konum belirtiyoruz.
     ToplamKarincaSayisi = round(sqrt(goruntuSatir*goruntuSutun));   %karakokunu aliyoruz ardindan degeri yuvarliyoruz.
     KarincaninKonumu = zeros(ToplamKarincaSayisi, 2); 

%% karincalarin satir ve sutun olarak konumlari belirleniyor.
     rand('state', sum(clock));
     geciciToplamKarincaSayisi = rand(ToplamKarincaSayisi, 2);
     KarincaninKonumu(:,1) = round(1 + (goruntuSatir-1) * geciciToplamKarincaSayisi(:,1));  
     KarincaninKonumu(:,2) = round(1 + (goruntuSutun-1) * geciciToplamKarincaSayisi(:,2));  
     arananYol = '8';

%% karincalarin hafiza/bellek uzunluklari tanimlanmaktadir. goruntunun
% satir sutun degerlerindeki boyutlara gore degerler atariz ve sonrasinda
% bir sonraki karincaninda bu bilgiye ulasabilmesi icin karincanin
% hafizasina degeri atariz.
      if goruntuSatir*goruntuSutun == 128*128    
        deger = 40;
      elseif goruntuSatir*goruntuSutun == 256*256
        deger = 30;
      elseif goruntuSatir*goruntuSutun == 512*512
        deger = 20;
      end               
      hafizaUzunlugu = round(rand(1).*(1.15*deger-0.85*deger)+0.85*deger);

%% karincanin bellegindeki satir ve sutunlara gore pozisyonlari kaydederiz.
     KarincaHafizasi = zeros(ToplamKarincaSayisi, hafizaUzunlugu);   

%% kaydettigimiz satir ve sutun bilgilerini goruntunun boyutlarina gore
% ayarlayarak karincalarin toplamda ne kadar adim atacaklarini belirliyoruz.
% 2 boyutlu matrislik verileri tek boyuta donusturuyoruz bununla birlikte 
% ne kadar tur atacaklarina dair iterasyon sayisini belirliyoruz. 
         ToplamAdimSayisi = str2num(get(handles.edit2, 'String')); 
      
%% karinca hafizasindaki bellek uzunluguna gore iterasyondaki sutun ve
% satir olarak pozisyonunu kaydediyor. 
% deltaferomon ile onceki karincanin hangi konumda olduguna dair satir
% sutun verilerini pozisyon olarak bir vektorde tutuyor.
      for IterasyonSayisi = 1:ToplamIterasyonSayisi   
        deltaFeromon = zeros(goruntuSatir, goruntuSutun);        
        for AdimSayisi = 1:ToplamAdimSayisi  
            secilenDeltaFeromon = zeros(goruntuSatir, goruntuSutun); 
            for SecilenKarinca = 1:ToplamKarincaSayisi   
                SecilenKarincaSatir = KarincaninKonumu(SecilenKarinca,1); 
                SecilenKarincaSutun = KarincaninKonumu(SecilenKarinca,2);
                % karincanin gecerli konumuna ait etrafindaki konumlari
                % degiskenlere aktariyor.
                if arananYol == '4'
                    satir = SecilenKarincaSatir;
                    sutun = SecilenKarincaSutun;
                    geciciYolAraliklari = [satir-1 sutun; satir sutun+1; 
                                            satir+1 sutun; satir sutun-1];
                elseif arananYol == '8'
                    satir = SecilenKarincaSatir;
                    sutun = SecilenKarincaSutun;
                    geciciYolAraliklari = [satir-1 sutun-1; satir-1 sutun; satir-1 sutun+1; 
                                             satir sutun-1; satir sutun+1; satir+1 sutun-1;
                                             satir+1 sutun; satir+1 sutun+1];
                end

                % goruntuye ait gecerli pozisyon araliklarini bulup gecici 
                % olarak bilgileri sakliyoruz.  
                geciciToplamKarincaSayisi = find(geciciYolAraliklari(:,1)>=1 & geciciYolAraliklari(:,1)<=goruntuSatir & geciciYolAraliklari(:,2)>=1 & geciciYolAraliklari(:,2)<=goruntuSutun);
                yolAraliklari = geciciYolAraliklari(geciciToplamKarincaSayisi, :);

                % gecici verilere ait konumlarin olasiliklari hesaplaniyor.    
                karincaninGecisOlasiligi_v = zeros(size(yolAraliklari,1),1);     
                karincaninGecisOlasiligi_p = zeros(size(yolAraliklari,1),1);     

                % yol araliklari karincanin hafizasinda olup olmadigina
                % dair sorgulamalari yapiyoruz. goruntuyu ve feromon
                % bilgilerini aktariyoruz.
                for i = 1:size(yolAraliklari,1)
                    geciciToplamKarincaSayisi = (yolAraliklari(i,1)-1)*goruntuSutun + yolAraliklari(i,2);
                    if length(find(KarincaHafizasi(SecilenKarinca,:)==geciciToplamKarincaSayisi))==0  % karincanin hafizasindaki konumlar degil ise
                        karincaninGecisOlasiligi_v(i) = goruntuKopyasi(yolAraliklari(i,1), yolAraliklari(i,2));
                        karincaninGecisOlasiligi_p(i) = feromon(yolAraliklari(i,1), yolAraliklari(i,2));
                    else % karincanin hafizasindaki konumlar ise
                        karincaninGecisOlasiligi_v(i) = 0;
                        karincaninGecisOlasiligi_p(i) = 0;                    
                    end                    
                end
                
%% etraftaki konumlar karincanin hafizasindaki verilerle ortusuyorsa
% tekrardan yeni kisa yol arayisina gitmek icin yeniden hesaplanmasi
% gerekiyor.
                if (sum(sum(karincaninGecisOlasiligi_v))==0) | (sum(sum(karincaninGecisOlasiligi_p))==0)                
                    for i = 1:size(yolAraliklari,1)
                        geciciToplamKarincaSayisi = (yolAraliklari(i,1)-1)*goruntuSutun + yolAraliklari(i,2);
                        karincaninGecisOlasiligi_v(i) = goruntuKopyasi(yolAraliklari(i,1), yolAraliklari(i,2));
                        karincaninGecisOlasiligi_p(i) = feromon(yolAraliklari(i,1), yolAraliklari(i,2));
                    end
                end                        

                % karinca algoritmasinin esas olasilik formulunu
                % kullanarak karincanin gecislerini hesapliyoruz.
                karincaGecisOlasiligi = (karincaninGecisOlasiligi_v.^alpha) .* (karincaninGecisOlasiligi_p.^beta) ./ (sum(sum((karincaninGecisOlasiligi_v.^alpha) .* (karincaninGecisOlasiligi_p.^beta))));       

                % karincanin sonraki adiminda rastgeleligi belirlemek icin
                % rand ile konumunu belirleriz.
                rand('state', sum(100*clock));     
                geciciToplamKarincaSayisi = find(cumsum(karincaGecisOlasiligi)>=rand(1), 1);

                SiradakiKarincaSatir = yolAraliklari(geciciToplamKarincaSayisi,1);   
                SiradakiKarincaSutun = yolAraliklari(geciciToplamKarincaSayisi,2); 

                if length(SiradakiKarincaSatir) == 0
                    SiradakiKarincaSatir = SecilenKarincaSatir;
                    SiradakiKarincaSutun = SecilenKarincaSutun;
                end

                KarincaninKonumu(SecilenKarinca,1) = SiradakiKarincaSatir;
                KarincaninKonumu(SecilenKarinca,2) = SiradakiKarincaSutun;

                secilenDeltaFeromon(KarincaninKonumu(SecilenKarinca,1), KarincaninKonumu(SecilenKarinca,2)) = 1;

                % karincanin hafizasina yeni konumlari atiyoruz. 
                if AdimSayisi <= hafizaUzunlugu
                    KarincaHafizasi(SecilenKarinca,AdimSayisi) = (KarincaninKonumu(SecilenKarinca,1)-1)*goruntuSutun + KarincaninKonumu(SecilenKarinca,2);
                elseif AdimSayisi > hafizaUzunlugu
                    KarincaHafizasi(SecilenKarinca,:) = circshift(KarincaHafizasi(SecilenKarinca,:),[0 -1]); % dairesel sekilde bir dizi olusturuyoruz.
                    KarincaHafizasi(SecilenKarinca,end) = (KarincaninKonumu(SecilenKarinca,1)-1)*goruntuSutun + KarincaninKonumu(SecilenKarinca,2);
                end
                
%% feromon olasilik fonksiyonuna simdiye kadar uretmis oldugumuz tum
% verileri yukleyerek artik her bir karincaya ait konumlari ilgili
% fonksiyona gonderiyoruz. 
                feromon = ((1-rho).*feromon + rho.*secilenDeltaFeromon.*goruntuKopyasi).*secilenDeltaFeromon + feromon.*(abs(1-secilenDeltaFeromon));
            end 
            
%% en son feromon'a ait veriyi gunceleyerek feromon sivi miktarini
% buharlasmalara karsi yeniden guncelliyoruz.
            deltaFeromon = (deltaFeromon + (secilenDeltaFeromon>0))>0;
            feromon = (1-phi).*feromon;   
        end % karincaya ait adim sayisini sonlandiriyoruz.
      end % iterasyonlari tamamliyoruz.
      
%% image'larin histogramlari hesaplanarak goruntu haline donusturulebilmesi
% icin Threshold degerleri hesaplatiyoruz.
    hesapla = histogramAyarlama(feromon);  
    set(handles.text6, 'String', 'Image islendi...');  
    imwrite(uint8(abs((feromon>=hesapla).*255-255)), gray(256), [dosyaAdi '_kkaKenarTespit_.bmp'], 'bmp');    
   
%% olusan image kaydedildikten sonra kullaniciya gui ile gosteriliyor.
        goster = handles.axes2;
        axes(goster);
        imshow([dosyaAdi '_kkaKenarTespit_.bmp']);        
           
         set(handles.text11, 'String', 'Toplam sure : ');
        set(handles.text7, 'String', toc); 
        set(handles.text12, 'String', 'Saniyedir.');
  toc
    
  
function gonder = histogramAyarlama(goruntu)

goruntu = goruntu(:);

%% goruntunun histogram ile birlikte yogunlugu hesapliyoruz T=mean(I)
[sayac, N]=hist(goruntu,256);
i=1;
kumulatif=cumsum(sayac);   % kumulatif toplam yapiyoruz.
T(i)=(sum(N.*sayac))/kumulatif(end);

% T'nin alt ve ust degerleri hesaplaniyor.
kumulatif2=cumsum(sayac(N<=T(i)));
MBT=sum(N(N<=T(i)).*sayac(N<=T(i)))/kumulatif2(end);

kumulatif3=cumsum(sayac(N>T(i)));
MAT=sum(N(N>T(i)).*sayac(N>T(i)))/kumulatif3(end);
i=i+1;
T(i)=(MAT+MBT)/2;

% T(i)~=T(i-1)
Threshold=T(i);
while abs(T(i)-T(i-1))>=1
    kumulatif2=cumsum(sayac(N<=T(i)));
    MBT=sum(N(N<=T(i)).*sayac(N<=T(i)))/kumulatif2(end);
    
    kumulatif3=cumsum(sayac(N>T(i)));
    MAT=sum(N(N>T(i)).*sayac(N>T(i)))/kumulatif3(end);
    
    i=i+1;
    T(i)=(MAT+MBT)/2; 
    Threshold=T(i);
end

gonder = Threshold;

function slider1_Callback(hObject, eventdata, handles)


function slider1_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function slider2_Callback(hObject, eventdata, handles)
global alpha 
alpha = get(hObject,'Value'); 
alpha = num2str(alpha); 
set(handles.edit3, 'String', alpha); 

function slider2_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

 
function slider3_Callback(hObject, eventdata, handles)
global beta  
beta = get(hObject,'Value'); 
set(handles.edit4, 'String', beta); 


function slider3_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function slider4_Callback(hObject, eventdata, handles)
global rho
rho = get(hObject,'Value'); 
set(handles.edit5, 'String', rho); 

function slider4_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function slider5_Callback(hObject, eventdata, handles)
global phi 
phi = get(hObject,'Value'); 
set(handles.edit6, 'String', phi); 


function slider5_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit3_Callback(hObject, eventdata, handles)

function edit3_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit4_Callback(hObject, eventdata, handles)

function edit4_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit5_Callback(hObject, eventdata, handles)

function edit5_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit6_Callback(hObject, eventdata, handles)

function edit6_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function figure1_CreateFcn(hObject, eventdata, handles)


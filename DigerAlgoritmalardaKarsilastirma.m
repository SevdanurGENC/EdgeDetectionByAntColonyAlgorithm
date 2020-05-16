
goruntu = imread('KameraliAdam.bmp');  
robertFiltresi = edge(goruntu,'roberts');
prewittFiltresi = edge(goruntu,'prewitt');
sobelFiltresi = edge(goruntu,'sobel');
LoGFiltresi = edge(goruntu,'log');   
cannygoruntu = edge(goruntu,'canny'); 

goruntuANT = imread('KameraliAdam.bmp_kkaKenarTespit_.bmp'); 

figure 
% subplot(241), imshow(goruntu), title('Orjinal Goruntu')
% subplot(242), imshow(1-robertFiltresi) , title('Robert Cross Filtresi')
% subplot(243), imshow(1-prewittFiltresi), title('Prewitt Filtresi')
% subplot(244), imshow(1-sobelFiltresi), title('Sobel Filtresi')
% subplot(245), imshow(1-LoGFiltresi), title('LoG Filtresi')
% subplot(246), imshow(1-cannygoruntu), title('Canny Filtresi')
% subplot(247), imshow(goruntuANT), title('ANT Filtresi')

subplot(251), imshow(goruntu), title('Orjinal Goruntu')
subplot(252), imshow(goruntuANT), title('ANT Filtresi')
%subplot(242), imshow(1-robertFiltresi) , title('Robert Cross Filtresi')
subplot(253), imshow(1-prewittFiltresi), title('Prewitt Filtresi')
subplot(254), imshow(1-sobelFiltresi), title('Sobel Filtresi')
%subplot(245), imshow(1-LoGFiltresi), title('LoG Filtresi')
subplot(255), imshow(1-cannygoruntu), title('Canny Filtresi')



clc;
clear;
close all;
imtool close all;
workspace;

image = '';
prompt = 'Enter an image filename: ';

while (isempty(image))
        image = input(prompt,'s');
end

prompt = 'Enter amount of characters: ';
amt = input(prompt);

%Resize the image
image = strcat(image,'.jpg');
carImage = imresize(imread(image),[640,NaN]);
figure(1);
imshow(carImage);
%Convert the image to gray_Scale
grayImage = rgb2gray(carImage);
[rows, cols] = size(grayImage);
dilateImage = grayImage;
%The for loop below performs the dialation of the image
for i = 1:rows
    for j = 2:cols-1
        temp = max(grayImage(i,j-1),grayImage(i,j));
        dilateImage(i,j) = max(temp,grayImage(i,j+1));
    end
end

image = dilateImage;
figure(2);
imshow(grayImage);
figure(3);
title('Dilated Image');
imshow(dilateImage);

figure(4);
imshow(image);

difference = 0;
sum = 0;
total_sum = 0;
difference = uint32(difference);

%% Horizontal Edge Detection
% The histogram represents column wise sum of the difference in the
% grayscale value of the neighbouring pixels values
disp('Processing Edges Horizontally...');
horz_max = 0;
maximum = 0;
for i = 2:cols
    sum = 0;
    for j = 2:rows
        if (image(j,i) > image(j-1,i))
            difference = uint32(image(j,i) - image(j-1,i));
        else
            difference = uint32(image(j-1,i) - image(j,i));
        end
        
        if (difference > 20)
            sum = sum + difference;
        end
    end
    horz_a(i) = sum;
    
    if (sum > maximum)
        horz_max = i;
        maximum = sum;
    end
    
    total_sum = total_sum + sum;
end
average = total_sum/cols;
figure(5)
subplot(3,1,1);
plot(horz_a);
title('Horizontal Edge Processing Histogram');
xlabel('Column Number ->');
ylabel('Difference ->');

disp('Passing Horizontal Histogram through LPF...');
sum = 0;
horz = horz_a;
for i = 21:(cols-21)
    sum = 0;
    for j = (i-20):(i+20)
        sum = sum + horz_a(j);
    end
    horz(i) = sum/41;
end
subplot(3,1,2);
plot(horz);
title('Histogram after passing through LPF');
xlabel('Column Number ->');
ylabel('Difference ->');

disp('Filter out Horizontal Histogram...');
for i = 1:cols
    if (horz(i) < average)
        horz(i) = 0;
        for j = 1:rows
            image(j,i) = 0;
        end
    end
end

subplot(3,1,3);
plot(horz);
title('Histogram after Filtering');
xlabel('Column Number ->');
ylabel('Difference ->');

difference = 0;
total_sum = 0;
difference = uint32(difference);

%% Vertical Edge Detection
% The histogram represents row wise sum of the difference in the
% grayscale value of the neighbouring pixels values
disp('Processing Edges Vertically...');
maximum = 0;
max_vert = 0;
for i = 2:rows
    sum = 0;
    for j = 2:cols
        if (image(i,j) > image(i,j-1))
            difference = uint32(image(i,j) - image(i,j-1));
        end
        if (image(i,j) <= image(i,j-1))
            difference = uint32(image(i,j-1) - image(i,j));
        end
        
        if (difference > 20)
            sum = sum + difference;
        end
    end
    
    vert_a(i) = sum;
    if (sum > maximum)
        vert_max = i;
        maximum = sum;
    end
    total_sum = total_sum + sum;
end
average = total_sum/rows;

figure(6);
subplot(3,1,1);
plot(vert_a);
title('Vertical Edge Processing Histogram...');
xlabel('Row Number ->');
ylabel('Difference ->');

%% Passing Vertical Histogram through LPF
disp('Passing Vertical Histogram through LPF...');
sum = 0;
vert = vert_a;

for i = 21:(rows-21);
    sum = 0;
    for j = (i-20):(i+20)
        sum = sum + vert_a(j);
    end
    vert(i) = sum/41;
end
subplot(3,1,2);
plot(vert);
title('Histogram after passing through LPF');
xlabel('Row Number ->');
ylabel('Difference ->');

disp('Filter out Vertical Histogram...');
for i = 1:rows;
    if (vert(i) < average)
        vert(i) = 0;
        for j = 1:cols
            image(i,j) = 0;
        end
    end
end
subplot(3,1,3);
plot(vert);
title('Histogram after Filtering');
xlabel('Row Number ->');
ylabel('Difference ->');

figure(7);
imshow(image);

j = 1;
for i = 2:cols-2
    if ((horz(i) ~= 0) && (horz(i-1) == 0) && (horz(i+1) == 0))
        column(j) = i;
        column(j+1) = i;
        j = j + 2;
    elseif(((horz(i) ~= 0) && (horz(i-1) == 0)) || ((horz(i) ~= 0) && (horz(i+1) == 0)))
        column(j) = i;
        j = j + 1;
    end
end

j = 1;
for i = 2:rows-2
    if ((vert(i) ~= 0) && (vert(i-1) == 0) && (vert(i+1) == 0))
        row(j) = i;
        row(j+1) = i;
        j = j + 2;
    elseif(((vert(i) ~= 0) && (vert(i-1) == 0)) || ((vert(i) ~= 0) && (vert(i+1) == 0)))
        row(j) = i;
        j = j + 1;
    end
end

[temp, col_size] = size(column);
if (mod(col_size,2))
    column(col_size+1) = cols;
end

[temp, row_size] = size(row);
if (mod(row_size,2))
    row(row_size+1) = rows;
end

%% This is just to display the extracted license plate, nothing special :)
for i = 1:2:row_size
    for j = 1:2:col_size
        if(~((horz_max >= column(j) && horz_max <= column(j+1)) && (vert_max >= row(i) && vert_max <= row(i+1))))
            for m = row(i):row(i+1)
                for n = column(j):column(j+1)
                    image(m,n) = 0;
                end
            end
        end
    end
end

%% PREPROCESSING FOR SEGMENTATION
figure, imshow(image);
BW = im2bw(image,0.5); %converts image to binary with threshold of .5
figure, imshow(BW);
BW = imcomplement(BW); %inverts the image
figure, imshow(BW);
final=bwareaopen(BW,100); % removes any image with pixel area bigger than 100
figure, imshow(BW);
props = regionprops(BW,'Image', 'BoundingBox'); %find image properties
                     %Image containes all the images broken apart
                     %Bound box contain the corrdinates and size of all the images
imshow(props)

                     
%% MY SEGMENTATION
NR=cat(1,props.BoundingBox); % Listing bound box info into a matrix
[M,N] = size(NR);
NR = NR(2:M,4)';  %changes it to an array of the images y-axis size
                %removes the first image bound box because it is info of the whole image
                %removes everything except the y-axis size

[amount,bin]=hist(NR); % Histogram of the y-dimension widths of all boxes.
ind=find(amount==amt); % Find the index of the bin which has 'amt'(amount of character enterd in the begining) common y-axis lenght.

if  ind - 1 == 0 
    binLow = 0; %if index is the first bin, for code resilience
elseif ind + 1 > length(amount) 
    binHigh = 100; %if index is the last bin, for code resilience
else
    binLow = (bin(ind)+bin(ind-1))/2; %Define the low end of bin with 'amt' common y-axis lenght.
    binHigh = (bin(ind)+bin(ind+1))/2; %Define the high end of bin with 'amt' common y-axis lenght.
end

index = find(NR >= binLow & NR <= binHigh) + 1; %finds the indexs of all the boxes witch are
                                               %within the bin the contains
                                               %contains 'amt' comon y-axis
                                               %lenghts
                                               %Accounts for the removal of
                                               %the first image above
                                               
%% TESTING SEGMENTATION
image = {props.Image}; 
figure
for i = 1:length(index)
    subplot(1,length(index),i), imshow(image{index(i)})
end

%% CHARACTER EXTRACTION

if isempty(index)
    fprintf('Unable to read license plate');
else
    load NewTemplates % Loads the templates of characters in the memory.
    for i = 1:length(index)
        check=imresize(image{index(i)},[42 24]); % Resize the input image so it can be compared with the template's images.
        comp=[ ];
        for n=1:length(NewTemplates)
            sem=corr2(NewTemplates{1,n},check); % Correlation the input image with every image in the template for best matching.
            comp= [comp sem]; % Record the value of correlation for each template's character.
        end
        vd=find(comp==max(comp)); % Find the index which correspond to the highest matched character.
    %*-*-*-*-*-*-*-*-*-*-*-*-*-
    % Accodrding to the index assign to 'letter'.
        % Alphabets listings.
        if vd==1 || vd==2
            letter='A';
        elseif vd==3 || vd==4
            letter='B';
        elseif vd==5
            letter='C';
        elseif vd==6 || vd==7
            letter='D';
        elseif vd==8
            letter='E';
        elseif vd==9
            letter='F';
        elseif vd==10
            letter='G';
        elseif vd==11
            letter='H';
        elseif vd==12
            letter='I';
        elseif vd==13
            letter='J';
        elseif vd==14
            letter='K';
        elseif vd==15
            letter='L';
        elseif vd==16
            letter='M';
        elseif vd==17
            letter='N';
        elseif vd==18 || vd==19
            letter='O';
        elseif vd==20 || vd==21
            letter='P';
        elseif vd==22 || vd==23
            letter='Q';
        elseif vd==24 || vd==25
            letter='R';
        elseif vd==26
            letter='S';
        elseif vd==27
            letter='T';
        elseif vd==28
            letter='U';
        elseif vd==29
            letter='V';
        elseif vd==30
            letter='W';
        elseif vd==31
            letter='X';
        elseif vd==32
            letter='Y';
        elseif vd==33
            letter='Z';
            %*-*-*-*-*
        % Numerals listings.
        elseif vd==34
            letter='1';
        elseif vd==35
            letter='2';
        elseif vd==36
            letter='3';
        elseif vd==37 || vd==38
            letter='4';
        elseif vd==39
            letter='5';
        elseif vd==40 || vd==41 || vd==42
            letter='6';
        elseif vd==43
            letter='7';
        elseif vd==44 || vd==45
            letter='8';
        elseif vd==46 || vd==47 || vd==48
            letter='9';
        else
            letter='0';
        end
        plate(i) = letter;
    end
end

fprintf('\nLicense Plate Number: \n')
disp(plate);



    

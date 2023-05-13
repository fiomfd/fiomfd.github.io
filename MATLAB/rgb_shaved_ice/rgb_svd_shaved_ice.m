% load image file
I=imread('shaved_ice.jpg');
RGB=imresize(I, 1/5);
[R,G,B] = imsplit(RGB);
% imshow(RGB)
% whos R, whos G, whos B

% uint8 -> float64
DR=double(R);
% R010
[UR010,SR010,VR010]=svds(DR,010);
DR010=UR010*SR010*VR010';
R010=uint8(DR010);
% R030
[UR030,SR030,VR030]=svds(DR,30);
DR030=UR030*SR030*VR030';
R030=uint8(DR030);
% R050
[UR050,SR050,VR050]=svds(DR,50);
DR050=UR050*SR050*VR050';
R050=uint8(DR050);
% R070
[UR070,SR070,VR070]=svds(DR,070);
DR070=UR070*SR070*VR070';
R070=uint8(DR070);
% R090
[UR090,SR090,VR090]=svds(DR,090);
DR090=UR090*SR090*VR090';
R090=uint8(DR090);

% uint8 -> float64
DG=double(G);
% G010
[UG010,SG010,VG010]=svds(DG,10);
DG010=UG010*SG010*VG010';
G010=uint8(DG010);
% G030
[UG030,SG030,VG030]=svds(DG,30);
DG030=UG030*SG030*VG030';
G030=uint8(DG030);
% G050
[UG050,SG050,VG050]=svds(DG,50);
DG050=UG050*SG050*VG050';
G050=uint8(DG050);
% G070
[UG070,SG070,VG070]=svds(DG,070);
DG070=UG070*SG070*VG070';
G070=uint8(DG070);
% G090
[UG090,SG090,VG090]=svds(DG,090);
DG090=UG090*SG090*VG090';
G090=uint8(DG090);

% uint8 -> float64
DB=double(B);
% B010
[UB010,SB010,VB010]=svds(DB,10);
DB010=UB010*SB010*VB010';
B010=uint8(DB010);
% B030
[UB030,SB030,VB030]=svds(DB,30);
DB030=UB030*SB030*VB030';
B030=uint8(DB030);
% B050
[UB050,SB050,VB050]=svds(DB,50);
DB050=UB050*SB050*VB050';
B050=uint8(DB050);
% B070
[UB070,SB070,VB070]=svds(DB,070);
DB070=UB070*SB070*VB070';
B070=uint8(DB070);
% B090
[UB090,SB090,VB090]=svds(DB,090);
DB090=UB090*SB090*VB090';
B090=uint8(DB090);

% R+G+B
% RGB010
RGB010(:,:,1)=(R010);
RGB010(:,:,2)=(G010);
RGB010(:,:,3)=(B010);
% RGB030
RGB030(:,:,1)=(R030);
RGB030(:,:,2)=(G030);
RGB030(:,:,3)=(B030);
% RGB050
RGB050(:,:,1)=(R050);
RGB050(:,:,2)=(G050);
RGB050(:,:,3)=(B050);
% RGB070
RGB070(:,:,1)=(R070);
RGB070(:,:,2)=(G070);
RGB070(:,:,3)=(B070);
% RGB090
RGB090(:,:,1)=(R090);
RGB090(:,:,2)=(G090);
RGB090(:,:,3)=(B090);

figure;
subplot(2,3,1);
imshow(RGB);
title('Original R, G and B are full rank 448 by 794 matrices respectively.');
subplot(2,3,2);
imshow(RGB010);
title('Rank 10 Approximation');
subplot(2,3,3);
imshow(RGB030);
title('Rank 30 Approximation');
subplot(2,3,4);
imshow(RGB050);
title('Rank 50 Approximation');
subplot(2,3,5);
imshow(RGB070);
title('Rank 70 Approximation');
subplot(2,3,6);
imshow(RGB090);
title('Rank 90 Approximation');

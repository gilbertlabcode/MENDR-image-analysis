function FinalImage = readTiffStack(FileTif)

% Auxilliary function for loading a multi-slice tiff file into MATLAB.
% Copied from http://www.matlabtips.com/how-to-load-tiff-stacks-fast-really-fast/

InfoImage=imfinfo(FileTif);
mImage=InfoImage(1).Width;
nImage=InfoImage(1).Height;
NumberImages=length(InfoImage);
FinalImage=zeros(nImage,mImage,NumberImages,'uint16');

TifLink = Tiff(FileTif, 'r');
for i=1:NumberImages
   TifLink.setDirectory(i);
   FinalImage(:,:,i)=TifLink.read();
end
TifLink.close();

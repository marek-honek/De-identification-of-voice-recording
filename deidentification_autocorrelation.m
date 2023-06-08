clc
clear
close all

% Název souboru
speaker = 'Zuzka'
filename = append(speaker, '_original.m4a');

% Načtení zvukového záznamu
[data, fs] = audioread(filename);

% Převedení záznamu na mono
data_mono = data(:,1)+data(:,2);

% Odříznutí nul ze začátku
firstNonzeroIndex = find(data_mono, 1, 'first');
data_mono = data_mono(firstNonzeroIndex:end);

% normovat peak2peak na rozkmit -1:1
maxSignal = max(abs(data_mono));
data_mono = data_mono./maxSignal;

% Hloubka LPC
LpcDepthLongTerm = 20;
LpcDepthSegment = 50;

% Délka úseku v počtu vzorků
segmentLength = 1024;

% Délka úseku v sekundách
segmentDuration = segmentLength/fs;

% Počet úseků
numSegments = floor(length(data_mono) / segmentLength);

% Imaginární jednotka
j = 1i;

% Coeficient indexu frekvencí
freqCoef = (0:fs/segmentLength:fs/2);
freqCoef(segmentLength/2+2:segmentLength) = fliplr(freqCoef(2:segmentLength/2));

% Rozsekání zvukového záznamu na úseky
segments = zeros(numSegments, segmentLength);

for i = 1:numSegments
    startIndex = (i-1) * segmentLength + 1;
    endIndex = i * segmentLength;
    segments(i,:) = data_mono(startIndex:endIndex);
end

disp('rozsekáno')

% Prealokace proměnných
segmentAutoccorrelation = zeros(numSegments, 2*segmentLength-1);
segmentFrequencies = zeros(numSegments, segmentLength);

% Spektra po úsecích
for i=1:numSegments
    segmentFrequencies(i,:) = fft(segments(i,:));
    segmentAutoccorrelation(i,:) = xcorr(segments(i,:));
end

disp('spektra a autokorelace segmentů')
    
% Průměrování autokorelační funkce
averageAutoccorrelation = mean(segmentAutoccorrelation,1);

% dlouhodobé LPC koeficienty
averageLpcCoef(:) = levinson(averageAutoccorrelation(segmentLength:segmentLength+ LpcDepthLongTerm+1),LpcDepthLongTerm);

disp('dlouhodobé LPC koeficienty')

% Prealokace
deidentifiedSegmentAutoccorrelation = zeros(numSegments,LpcDepthSegment+1);
Ra = zeros(1,LpcDepthLongTerm+1);

% Rovnice 15
for m = 0:LpcDepthLongTerm
    for i = 0:LpcDepthLongTerm-m
        Ra(m+1) = Ra(m+1) + averageLpcCoef(i+1) * averageLpcCoef(i+m+1);
    end
end

% Rovnice 14
for s = 1:numSegments
    for k  = 0:LpcDepthSegment+1
        sum = 0;
        for m = 1:LpcDepthLongTerm
            sum = sum + Ra(m+1) * (segmentAutoccorrelation(s,segmentLength+abs(k-m)) + segmentAutoccorrelation(s,segmentLength+abs(k+m)));
        end
        deidentifiedSegmentAutoccorrelation(s, k+1) = (Ra(1)*segmentAutoccorrelation(s,segmentLength)+sum);
    end
end

disp('deidentifikované autokorelační funkce')

% Prealokace
deidentifiedSegmentLpcCoef = zeros(numSegments,LpcDepthSegment+1);
deidentifiedSegmentLpcSpectrum = zeros(numSegments,segmentLength);

% LPC coef po segmentech
for i = 1:numSegments
    deidentifiedSegmentLpcCoef(i,:) = levinson(deidentifiedSegmentAutoccorrelation(i,1:LpcDepthSegment+2),LpcDepthSegment);
end

% LPC spectrum po segmentech
for i = 1:numSegments
    for f = 1:segmentLength
        sum = 0;    
        for sumIndex = 1:LpcDepthSegment
            sum = sum + deidentifiedSegmentLpcCoef(i,sumIndex+1) * exp(-j * sumIndex * 2 * pi * freqCoef(f) / fs);
        end
        deidentifiedSegmentLpcSpectrum(i,f) = 1/abs(1+sum);
    end
end

disp('deidentifikovaná LPC spektra')

% Prealokace
deidentifiedSegmentFrequencies = zeros(numSegments, segmentLength);

% Součet reálné složky LPC spektra a imaginární složky FFT spektra 
for i=1:numSegments
    deidentifiedSegmentFrequencies(i,:) = deidentifiedSegmentLpcSpectrum(i,:) .* cos(angle(segmentFrequencies(i,:))) + j .* deidentifiedSegmentLpcSpectrum(i,:) .* sin(angle(segmentFrequencies(i,:)));
end

% Inicializace proměnné pro upravené segmenty
deidentifiedSegment = zeros(numSegments, segmentLength);

% Deidentifikované segmenty
for i=1:numSegments
    deidentifiedSegment(i,:) = ifft(deidentifiedSegmentFrequencies(i,:));
end

% Prealokace
deidentifiedSegmentFrequencies = zeros(numSegments, segmentLength);

% Součet reálné složky LPC spektra a imaginární složky FFT spektra 
for i=1:numSegments
deidentifiedSegmentFrequencies(i,:) = deidentifiedSegmentLpcSpectrum(i,:).*cos(angle(segmentFrequencies(i,:)))+j.*deidentifiedSegmentLpcSpectrum(i,:).*sin(angle(segmentFrequencies(i,:)));
end

% Inicializace proměnné pro upravené segmenty
deidentifiedSegmentFromLPC = zeros(numSegments, segmentLength);

% Deidentifikované segmenty
for i=1:numSegments
    deidentifiedSegmentFromLPC(i,:) = real(ifft(deidentifiedSegmentFrequencies(i,:)));
end

% Srovnání středních hodnot
meanValue = mean(deidentifiedSegmentFromLPC,2);
deidentifiedSegmentFromLPC = deidentifiedSegmentFromLPC-meanValue;

% Srovnání hlasitosti segmetů
meanSegmentSource = mean(abs(segments),2);
meanSegmentOutput = mean(abs(deidentifiedSegmentFromLPC),2);
deidentifiedSegmentFromLPC = deidentifiedSegmentFromLPC.*(meanSegmentSource./meanSegmentOutput);

disp('deidentifikované segmenty')

% Vyhlazení okrajů segmentů
window = ones(1,segmentLength);
window(1) = 0.001;
window(2) = 0.002;
window(3) = 0.005;
window(4) = 0.01;
window(5) = 0.02;
window(6) = 0.05;
window(7) = 0.1;
window(8) = 0.2;
window(9) = 0.5;
window(segmentLength) = 0.001;
window(segmentLength-1) = 0.002;
window(segmentLength-2) = 0.005;
window(segmentLength-3) = 0.01;
window(segmentLength-4) = 0.02;
window(segmentLength-5) = 0.05;
window(segmentLength-6) = 0.1;
window(segmentLength-7) = 0.2;
window(segmentLength-8) = 0.5;
deidentifiedSegmentFromLPC = deidentifiedSegmentFromLPC.*window;

disp('vyhlazené přechody mezi segmenty')

% Inicializace proměnné pro sestavení časového průběhu
reconstructedSignal = zeros(1, segmentLength * numSegments);

% Přidání rekonstruovaného úseku do celkového časového průběhu
for i = 1:numSegments
    startIndex = (i - 1) * segmentLength + 1;
    endIndex = i * segmentLength;
    reconstructedSignal(startIndex:endIndex) = deidentifiedSegmentFromLPC(i,:);
end

% Filtrování nízkých frekvencí
wholeSignalSpectrum = fft(reconstructedSignal);
BPfiter = zeros(1,length(reconstructedSignal));
for i=round(length(reconstructedSignal)/240):round(length(reconstructedSignal)/12)
    BPfiter(i)=1;
end

for i=round(11*length(reconstructedSignal)/12):round(239*length(reconstructedSignal)/240)
    BPfiter(i)=1;
end
wholeSignalSpectrum=wholeSignalSpectrum.*BPfiter;
filteredReconstructedSignal = real(ifft(wholeSignalSpectrum));

% Normování peak2peak na rozkmit -1:1
maxSignal = max(abs(filteredReconstructedSignal));
filteredReconstructedSignal = filteredReconstructedSignal./maxSignal;

disp('hotovo')

figure
plot(data_mono(:))
hold on
plot(filteredReconstructedSignal(:))

sound(filteredReconstructedSignal, fs)

filename = append(speaker, '_' , num2str(LpcDepthLongTerm),'_', num2str(LpcDepthSegment), '.m4a');
audiowrite(filename, filteredReconstructedSignal, fs)

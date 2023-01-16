function raw = retro_gating_fill(raw, phaseNum)
%
%   -- Lamy Jerome -- 2019/07, last update: 2022/02
%   jimolrame@gmail.com

display('Starting Filling Retro Gating...')
tic
% phaseNum = find(strcmp(header,'006:sMDH.sLC.ushPhase')) ;
s = size(raw) ;

lDim = 1 : numel(s) ;
lDimParse = lDim ;
lDimParse(phaseNum) = [] ;

raw = reshape( permute( raw, [ phaseNum, lDimParse ] ), s(phaseNum), numel(raw)/s(phaseNum) ) ;

sR = size(raw) ;

id =  raw(end,:)==0 & raw(1,:)~=0  ;
l = sum(id) ;

%-- If no end missing data
if l==0
    raw = permute( reshape( raw, [ s(phaseNum), s(lDimParse) ]), [lDim(2:phaseNum) 1 lDim(phaseNum+1:end) ] )  ;
    display('No retrogating missing data')
    return ;
end

% only interpolate from xx phase :
g = .65.* (1-sum(min(abs(raw(:,raw(1,:)~=0)'))==0)/size(raw,1)) ; % default 65% but some of the patients have really bad rythm

if g==0
    raw_keep = [] ;
else
    raw_keep = raw(1:floor(g*size(raw,1)),:);
    raw = raw(floor(g*size(raw,1))+1:end,:) ;
end

% -- id = elimination *
nRaw = raw ;
for k = 2 : size(raw,1)
    idInterp = (raw(1,:)~=0) & (sum(abs(nRaw(k:end,:)),1)==0) ;
   
    sId = numel(find(idInterp)) ;
    if sId == 0
       continue
    end
    sMakima = interp1(...
        1:k-1,...
        [ real(nRaw(1:k-1,idInterp)) , imag(nRaw(1:k-1,idInterp)) ],...
        1: (k-2)/(size(raw,1)-1)/100 :k-1,...
        'pchip') ;
    sMakima = sMakima(1:100:end,:) ;

    nRaw(:, idInterp) = sMakima(:,1:sId) + 1i * sMakima(:,sId+1:end) ;
    display([num2str(sId/size(raw,2)*100) '% columns were missing ' num2str((size(raw,1)-k+1)/sR(1)*100) '% information'])
end

raw = nRaw ;



%%
raw = [raw_keep;raw] ;


raw = permute( reshape( raw, [ s(phaseNum), s(lDimParse) ]), [lDim(2:phaseNum) 1 lDim(phaseNum+1:end) ] )  ;

display(['Filling Retro Gating done in = ' num2str(toc) ' s'])

end
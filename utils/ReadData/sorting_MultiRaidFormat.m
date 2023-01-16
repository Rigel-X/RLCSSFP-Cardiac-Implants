function [ raw, header, mr_times, mr_times_header ]  = sorting_MultiRaidFormat( MultiRAID )

%JL

% rawdata         = MultiRAID.rawdata;
% loopcounters    = MultiRAID.loopcounters;
% sMDH            = MultiRAID.sMDH;
% MrProt          = MultiRAID.MrProt;
% Yaps            = MultiRAID.Yaps;
% lc_names        = MultiRAID.lc_names;
% sCHA            = MultiRAID.sCHA;
% RaidInfo        = MultiRAID.RaidInfo;

protname = deblank(MultiRAID.RaidInfo.protname') ;
display(protname)

%lc_names ?
% 01 ->lines
% 02 -> Averaging ?? It is supposed to be
% 03 -> Slice = position (ex sax mutlitple ?)
% 04 -> Partition ???????? What is it ?
% 05 -> Echo ???????? What is it ? acquiered with different TE ?
% 06 -> Phase  = time phase ?
% 07 -> Repetition = Num if averaged ?
% 08 -> ushSet phi0 ,1 ,2 , 3 for the differences
% 09 -> ushSeg ???????? What is it ?
% 10 -> ushIda ???????? What is it ?
% 11 -> ushIdb ???????? What is it ?
% 12 -> ushIdc ???????? What is it ?
% 13 -> ushIdd ???????? What is it ?
% 14 -> ushIde ???????? What is it ?
% 15 -> ushChannel Coils number ?


%% HOW TO SORT DATA
%-- From this part we can deduce the lines encoded multiple times to choose
%from 019:sMDH.aulEvalInfoMask(1) and 020:sMDH.aulEvalInfoMask(2)
% course2follow = MultiRAID.lc_names' ;
% for k = 1 : length(MultiRAID.lc_names(1:64))
%     course2follow{2,k} = unique(MultiRAID.loopcounters(:,k)) ;
%     course2follow{3,k} = min(unique(MultiRAID.loopcounters(:,k))) ;
%     course2follow{4,k} = length(unique(MultiRAID.loopcounters(:,k))) ;
%     course2follow{5,k} = max(unique(MultiRAID.loopcounters(:,k))) ;
%
%     a = zeros(1,length( course2follow{2,k} )) ;
%     b = zeros(1,length( course2follow{2,k} )) ;
%     for ka = 1 : length( course2follow{2,k} )
%         a(ka) = numel( find( MultiRAID.loopcounters(:,k)==course2follow{2,k}(ka) ) ) ;
%         b(ka) = numel( find( MultiRAID.loopcounters(:,k)==course2follow{2,k}(ka) ) )/...
%             size(MultiRAID.loopcounters,1) ;
%     end
%     course2follow{7,k} = a ;
%
%     course2follow{8,k} = min(a) / max(a) ;
%
% end
% clearvars -except course2follow MultiRAID

%%
if ~isa(MultiRAID.rawdata, 'single')
    display('Attention raw data converted to complex single for memory management')
    MultiRAID.rawdata = single(MultiRAID.rawdata) ;
end

display('Initial matrix single complex')
raw = complex(zeros(...
    MultiRAID.sMDH.ushSamplesInScan,...
    max(MultiRAID.loopcounters(:,1)),...
    max(MultiRAID.loopcounters(:,2)),...
    max(MultiRAID.loopcounters(:,3)),...
    max(MultiRAID.loopcounters(:,4)),...
    max(MultiRAID.loopcounters(:,5)),...
    max(MultiRAID.loopcounters(:,6)),...
    max(MultiRAID.loopcounters(:,7)),...
    max(MultiRAID.loopcounters(:,8)),...
    max(MultiRAID.loopcounters(:,10)),... % 9 is seg, not needed
    max(MultiRAID.loopcounters(:,11)),...
    max(MultiRAID.loopcounters(:,12)),...
    max(MultiRAID.loopcounters(:,13)),...
    max(MultiRAID.loopcounters(:,14)),...
    max(MultiRAID.loopcounters(:,15)),'single')) ;

mr_times = zeros(...
    3,...
    max(MultiRAID.loopcounters(:,1)),...
    max(MultiRAID.loopcounters(:,2)),...
    max(MultiRAID.loopcounters(:,3)),...
    max(MultiRAID.loopcounters(:,4)),...
    max(MultiRAID.loopcounters(:,5)),...
    max(MultiRAID.loopcounters(:,6)),...
    max(MultiRAID.loopcounters(:,7)),...
    max(MultiRAID.loopcounters(:,8)),...
    max(MultiRAID.loopcounters(:,10)),... % 9 is seg, not needed
    max(MultiRAID.loopcounters(:,11)),...
    max(MultiRAID.loopcounters(:,12)),...
    max(MultiRAID.loopcounters(:,13)),...
    max(MultiRAID.loopcounters(:,14)),...
    max(MultiRAID.loopcounters(:,15))) ;
mr_times_header = MultiRAID.lc_names(30:32) ;

%% Missing data estimation
s = size(raw) ;
n = numel(raw) ;
[ avoidMultiple ] = missing_info(MultiRAID, s, n) ;

%% Filling matrix
display('Filling raw matrix')
tic
liste_k = 1 : size(MultiRAID.loopcounters,1) ;
liste_k(avoidMultiple) = [] ;
for k = liste_k
    raw(...
        :,...
        MultiRAID.loopcounters(k,1),...
        MultiRAID.loopcounters(k,2),...
        MultiRAID.loopcounters(k,3),...
        MultiRAID.loopcounters(k,4),...
        MultiRAID.loopcounters(k,5),...
        MultiRAID.loopcounters(k,6),...
        MultiRAID.loopcounters(k,7),...
        MultiRAID.loopcounters(k,8),...
        MultiRAID.loopcounters(k,10),... % 9 is seg, not needed
        MultiRAID.loopcounters(k,11),...
        MultiRAID.loopcounters(k,12),...
        MultiRAID.loopcounters(k,13),...
        MultiRAID.loopcounters(k,14),...
        MultiRAID.loopcounters(k,15)) = MultiRAID.rawdata(k,:) ;
    
    mr_times(...
        :,...
        MultiRAID.loopcounters(k,1),...
        MultiRAID.loopcounters(k,2),...
        MultiRAID.loopcounters(k,3),...
        MultiRAID.loopcounters(k,4),...
        MultiRAID.loopcounters(k,5),...
        MultiRAID.loopcounters(k,6),...
        MultiRAID.loopcounters(k,7),...
        MultiRAID.loopcounters(k,8),...
        MultiRAID.loopcounters(k,10),... % 9 is seg, not needed
        MultiRAID.loopcounters(k,11),...
        MultiRAID.loopcounters(k,12),...
        MultiRAID.loopcounters(k,13),...
        MultiRAID.loopcounters(k,14),...
        MultiRAID.loopcounters(k,15)) = MultiRAID.loopcounters(k,30:32) ;
    
end
display(['Sorting data time = ' num2str(toc) ' s'])

%% Squeeze + Header
header = [ {'000:sMDH.ushSamplesInScan'} ; MultiRAID.lc_names( [1:8,10:15] ) ] ;
s = size(raw) ;
header(s==1) = [] ;
raw = squeeze(raw) ;
mr_times = squeeze(mr_times) ;

end

function [ avoidMultiple ] = missing_info(MultiRAID, s, n)

display('Missing and double data evaluation')
tic
%-- evaluation of present data :
id = sub2ind(s(2:end) ,...
    MultiRAID.loopcounters(:,1),...
    MultiRAID.loopcounters(:,2),...
    MultiRAID.loopcounters(:,3),...
    MultiRAID.loopcounters(:,4),...
    MultiRAID.loopcounters(:,5),...
    MultiRAID.loopcounters(:,6),...
    MultiRAID.loopcounters(:,7),...
    MultiRAID.loopcounters(:,8),...
    MultiRAID.loopcounters(:,10),... % 9 is seg, not needed
    MultiRAID.loopcounters(:,11),...
    MultiRAID.loopcounters(:,12),...
    MultiRAID.loopcounters(:,13),...
    MultiRAID.loopcounters(:,14),...
    MultiRAID.loopcounters(:,15)) ;

%-- Num of not encoded
id_full = (1 : n/MultiRAID.sMDH.ushSamplesInScan)' ;
nZero = setdiff(id_full,id) ;

%-- Num of multiple encoded
[~,I,J] = unique(id, 'stable') ;

ixDupRows = setdiff((1:length(id))', I)' ;

[ id_dubs, ssId ] = sort(id([ J(ixDupRows); ixDupRows' ])) ;
lUnique = unique(id_dubs) ;
nMultiple = zeros(1,numel(lUnique+1)) ;
for k = 1 : numel(lUnique)
    nMultiple( numel(find(id_dubs==lUnique(k))) ) = nMultiple( numel(find(id_dubs==lUnique(k))) ) + 1 ;
end
if ~isempty(nMultiple)
    nMultiple(1) = [] ;
end
nMultiple(find(nMultiple,1,'last')+1:end) = [] ;
nMultiple = max(nMultiple) ;

%-- display
display( '#-------------------------------------------------------------#' )
display([ 'Total number of required lines : ' num2str(n/MultiRAID.sMDH.ushSamplesInScan)])
display([ 'Total number of encoded lines : ' num2str(numel(id))])
display([ 'Total number of unique lines : ' num2str(numel(unique(id))) ])
display([ 'Lines encoded 0 time : ' num2str( numel(nZero) ) ])
display([ 'Lines encoded 1 time : ' num2str(numel(unique(id)) - numel(id_dubs)) ])
display([ 'Lines encoded multiple times (2,3,...) : ' num2str(nMultiple') ' [ ]'])
display( '#-------------------------------------------------------------#' )

display(['Identifying missing or doubleb data time = ' num2str(toc)])

display('Choosing lines to keep from multiple encoding with EvalInfoMask(2)')
indexDubsId = [ J(ixDupRows); ixDupRows' ] ;

iSel = 20 ; % 020:sMDH.aulEvalInfoMask(2)
mDouble = MultiRAID.loopcounters(indexDubsId, :) ;

lSelId = unique( MultiRAID.loopcounters(:,iSel) ) ;
nSelId = zeros(size(lSelId)) ;
for k = 1 : numel(lSelId)
    nSelId(k) = numel(find( MultiRAID.loopcounters(:,iSel) == lSelId(k) )) ;
end

nSelDubs = zeros(size(lSelId)) ;
for k = 1 : numel(lSelId)
    nSelDubs(k) = numel(find( mDouble == lSelId(k) )) ;
end

[ lSel, idSelId, idSelDubs ] = intersect(lSelId, lSelId) ;
lSel = lSel( nSelId( idSelId ) == nSelDubs( idSelDubs ) ) ;
id = zeros(size(mDouble,1),1) ;
for k = 1 : length(lSel)
    id = id| (mDouble(:,iSel) == lSel(k)) ;
end
avoidMultiple = indexDubsId( id ) ;

if ~isempty(avoidMultiple)
    display('Discriminant found for doubled lines')
end


end

function fitswrite_bintable(data_in, info, filename)
%FITSWRITE saves a Matlab matrix as a FITS image
%
%Usage: fitswrite(data, filename)
%
%Known deficiencies
%
%1. Doesn't support arrays of imaginary numbers.
%2. Only handles simple 2 dimensional arrays of data.
%
%Author: R. G. Abraham, Institute of Astronomy, Cambridge University
%        abraham@ast.cam.ac.uk

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Modified 2012/08/01
% To write Binary Tables
%Following:
% www.cv.nrao.edu/fits/documents/standards/bintable_aa.ps
% (JET)

%e.g: >>data=fitsread(filename, 'BinTable');
%     >>info=fitsinfo(filename);
%     >>fitswrite_bintable(data, info, 'test.fits')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% try to figure out if we need to swap bytes. This is
% imperfect as I don't know the endian-ness of each
% architecture, so I'm only looking for ones I know for 
% sure are big-endian.
friend = computer;
if strmatch(friend,'PCWIN')
   bswap = 'b';
elseif strmatch(friend,'LNX86')
   bswap = 'b';   
elseif strmatch(friend,'ALPHA')
   bswap = 'b';
elseif strmatch(friend,'GLNXA64')
   bswap = 'b';  
else
   bswap = 'b';
end
%I made these all 'b' (JET)


%get hdr, xhdr
hdr=info.PrimaryData.Keywords;
infobin=info.BinaryTable;

%make hdr
header_cards=[];
for i=1:size(hdr,1)
  card=make_card(hdr{i,1}, hdr{i,2});
  header_cards=[header_cards; card];
end  

header_record = make_header_record(header_cards);
fid=fopen(filename,'w');
fwrite(fid,header_record','char');


%make xhdr

for k=1:length(infobin)
xhdr=infobin(k).Keywords;  
header_cards_bin=[];
for i=1:size(xhdr,1)
  card=make_card(xhdr{i,1}, xhdr{i,2});
  header_cards_bin=[header_cards_bin; card];
end  


header_record_bin = make_header_record(header_cards_bin);
fwrite(fid,header_record_bin','char');


data=data_in(:,:,k);


%----
if(iscell(data))
  dim=size(data,2);
  data=cell2mat(data);
else
  dim=size(data,2);
end




%find size of data
[nrow,ncol]=size(data);


%find out what type of data each dimension of data is
for j=1:dim
  
 %find what type of data this col is
 val=0;loc=1;
 while val == 0 && loc < length(xhdr)
   val=val+strcmp(xhdr{loc,1}, ['TFORM' num2str(j)]);
   loc=loc+1;
 end
 if loc==length(xhdr) %we didn't find tform specified
   error('Could not determine type of data to write')
 else 
   tform{j}=xhdr{loc-1,2};
   indiv_tform=tform{j};
 end
 
 %a few types it could be
 
 if(~isempty(strfind(tform{j}, 'E')))
   convert_tform{j}='float32';
   if(strfind(tform{j},'E') == 1) %in case the 1 is left off
     num(j)=1;
   else
     num(j)=str2num(indiv_tform(1:strfind(tform{j}, 'E')-1));
   end
 end

 if(~isempty(strfind(tform{j}, 'J')))
   convert_tform{j}='int32';
   if(strfind(tform{j},'J') == 1) %in case the 1 is left off
     num(j)=1;
   else
     num(j)=str2num(indiv_tform(1:strfind(tform{j}, 'J')-1));
   end
 end
 
 if(~isempty(strfind(tform{j}, 'D')))
   convert_tform{j}='float64';
   if(strfind(tform{j},'D') == 1) %in case the 1 is left off
     num(j)=1;
   else
     num(j)=str2num(indiv_tform(1:strfind(tform{j}, 'D')-1));
   end
 end
 
 if(~exist('convert_tform', 'var'))
   error('Could not determine type of data to write')
 end

end

% check if all the columns are of the same type. If
% yes, this will lead to a dramatic speed up in writing
% them to disk:
type_same = false;
if length(convert_tform)==1
  type_same=true;
else
  type_same = true;
  for ii=2:length(convert_tform)
    if ~strcmp(convert_tform(ii-1),convert_tform(ii))
      type_same=false;
    end
  end
end

% Put data in the right form
if size(data,1)==nrow & size(data,2)==ncol & type_same
  data=data';
end

% write the data, if all columns have the same type,
% avoid the loop and write all at once:
if type_same
  fwrite(fid,data, convert_tform{1}, bswap)
else
  % fits row/col is weird, so you have to write in a loop..
  for i=1:nrow
    st=1;en=st+num(1)-1;
    for j=1:dim
      fwrite(fid,data(i,st:en), convert_tform{j}, bswap);
      st=st+num(j);
      en=en+num(j);
    end
  end
end

%pad to integer * 2880
sz=ncol*nrow*4;
addn=2880-rem(sz,2880);
pad=ones(addn/4,1);
fwrite(fid,pad, 'float', bswap);

end  %multple bintables


fclose(fid);

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function card=make_card(keyword,value)
%MAKE_CARD turns a set of strings into a valid FITS card

%Make keyword field 8 bytes long
lk=length(keyword);
if (lk > 8) & (nargin>1)
	error('Keyword must be less than or equal to 8 characters!')
elseif (lk < 8 )
	keyword=[keyword,setstr(ones(1,8-lk)*32)];
end;

%Deal with both straight keyword and keyword/value pair
if (nargin==1 | strcmp(keyword,'END     '))
	%Keyword without a value
	card=keyword;	
else
	%Key/value pair has an equal sign and space at bytes 9 and 10
	card=[keyword,'= '];

	%Now output the value. The FITS standard wants things to start 
	%in different columns depending on what type of data the
	%value holds, according to the following rules:
	%
	%  Logical: T or F in column 30
	%
	%  Character string: A beginning quote in column 11 and an
	%  ending quote between columns 20 and 80.
	%
	%  Real part of an integer or floating point number: right 
	%  justified, ending in column 30.
	%
	%  Imaginary part: right justified, ending in
	%  column 50, starting after column 30 (NB. I won't bother 
	%  to support an imaginary part in this M-file, and will 
	%  let some radio astronomer who needs it add it if they want).

	if isstr(value)
  	    %Test for logical. If logical it goes in column 30 
		if (length(value)==1) & (strmatch(upper(value),'T') | strmatch(upper(value),'F'))
 			card=[card,setstr(ones(1,19)*32),value];	
		else	
			%Value must be a character string. Pad if less than 8
			%characters long.
			lv=length(value);
		    if (lv > 70)
		   error('Value must be less than 70 characters long!')
		    elseif (lv < 10 )
 	 	   value=[value,setstr(ones(1,8-lv)*32)];
		    end;
			card=[card,'''',value,''''];
		end;	
	else
		%Value must be a number. Convert to a string. Maximum
		%precision is set to 10 digits
		value=num2str(value,10);
		lv=length(value);
	
		%Write this out so it will end on column 30
		card=[card,setstr(ones(1,20-lv)*32),value];	
	end;
end;

%Now pad the output to make it exactly 80 bytes long
card=[card,setstr(ones(1,80-length(card))*32)];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function hrec=make_header_record(card_matrix)

[nrow,ncol] = size(card_matrix);
n_blanks = 36 - rem(nrow,36);
blank_line = setstr(ones(1,80)*32);
hrec = [card_matrix; repmat(blank_line,n_blanks,1)];
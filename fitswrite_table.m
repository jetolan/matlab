function fitswrite_table(data, hdr, xhdr, filename)
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
%Modified 2012/08/29
% To write Ascii Tables
%Following:
% www.cv.nrao.edu/fits/documents/standards/bintable_aa.ps
% (JET)

%e.g: >>data=fitsread(filename, 'table');
%     >>info=fitsinfo(filename);
%     >>hdr=info.PrimaryData.Keywords; xhdr=info.AsciiTable.Keywords;
%     >>fitswrite_table(data, hdr, xhdr, 'test.fits')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dim=size(data,2);

%convert to mat if cell
if(iscell(data))
  data=cell2mat(data);
end

%find size of data
[nrow,ncol]=size(data);

%make hdr
header_cards=[];
for i=1:size(hdr,1)
  card=make_card(hdr{i,1}, hdr{i,2});
  header_cards=[header_cards; card];
end  

%make xhdr
header_cards_bin=[];
for i=1:size(xhdr,1)
  card=make_card(xhdr{i,1}, xhdr{i,2});
  header_cards_bin=[header_cards_bin; card];
end  
    
    
header_record = make_header_record(header_cards);
header_record_bin = make_header_record(header_cards_bin);
%[ncards,dummy]=size(header_cards);
%fprintf(header_record(1,:));

fid=fopen(filename,'w');
fwrite(fid,header_record','char');
fwrite(fid,header_record_bin','char');





%Now write the table
%This is very specific to CAMB spec files at the moment.
%Other formats will need work

for i=1:size(data,1)
for j=1:size(data,2)
  if(j~=size(data,2))
    if(data(i,j)>=0)
    fprintf(fid,' %.7E  ', data(i,j));
    else
    fprintf(fid,'%.7E  ', data(i,j));
    end
  else
    if(data(i,j)>=0)
    fprintf(fid,' %.7E ', data(i,j));
    else
    fprintf(fid,'%.7E ', data(i,j));
    end
  end  
end
end


%pad to integer * 2880
dirf=dir(filename);
sz=dirf.bytes;
addn=2880-rem(sz,2880);
for i=1:addn
fprintf(fid, ' ');
end
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
if (nargin==1)
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
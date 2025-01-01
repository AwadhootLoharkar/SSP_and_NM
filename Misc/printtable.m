function printtable(a,varargin)
    
%% ROW-WISE PRINTING OF A NUMERICAL mXn ARRAY IN SCIENTIFIC NOTATION
%  thereby defining the precision as the number of significant digits.

%  NB: The concept of precision is not applicable to integer data.
    
%% INPUTS

%  a(1:m,1:n) = mXn array of real numbers to be printed row-wise
    
%% Examples of use 
    
%  Printing the row-vectors x,y and z, having the same length m,
%  using the default precision (6 significant digits):
%
%  >> printtable([x' y' z'])

%  Printing the row-vectors u,v,w,x,y and z (all having the same length)
%  stipulating the precisions:
%  u: 4 significant digits;
%  v: 3 significant digits;
%  w to z: 5 significant digits.
%
%  >> printtable([u',v',w',x',y',z'],'Precision',[4 3 5])

%  Printing the row-vectors u,v,w,x,y and z (all having the same length)
%  stipulating that w and y be printed in integer format
%  and that the other data be printed in normalized scientific notations:
%
%  >> printtable([u',v',w',x',y',z'],'Integer',[3 5],'Normalized',true)

%  Inserting line numbers requires the option 'Integer':
%    
%  >> printtable([ (1:m)' x' y' z'],'Integer',1)

%  Inserting line numbers starting at some integer value s :
%    
%  >> printtable([ (s:s+m)' x' y' z'],'Integer',1)

%% Recognized options in varargin 
% (uppercases for readability are optional): 

% if varargin{k} = 'Precision', 
% then varargin{k+1} = array p = [p(1) p(2) p(3) ...]
%                    = number of significant digits in the mantissae, i.e.:
%                      p(i) is the precision of a(:,i).
%                      If lp=length(p) < n=columns(a), 
%                      then p(lp) is the precision of a(:,lp:n).
%                      Default: p=6 for all columns of non integer data.
  
% if varargin{k} = 'Integer', 
% then varargin{k+1} = [i j k ...] = indices of columns 
%                      to be displayed as integers,
%                      thereby not applying the concept of precision
    
% if varargin{k} = 'Normalized', 
% then varargin{k+1} = logical (Boolean).
%                      If true, use normalized scientific notation 
    
% if varargin{k} = 'StartOfLine', 
% then varargin{k+1} = character string = separator before first field
    
% if varargin{k} = 'Separator', 
% then varargin{k+1} = character string = separator between fields

% if varargin{k} = 'EndOfLine', 
% then varargin{k+1} = character string = separator after last field
    
% if varargin{k} = 'European', 
% then varargin{k+1} = logical (Boolean). 
%                      If true, use European radix (comma instead
%                      of dot) in the output of non-integer numbers.
    
% if varargin{k} = 'CSV', 
% then varargin{k+1} = logical (Boolean). If true, use CSV table format.   

% if varargin{k} = 'HTML', 
% then varargin{k+1} = logical (Boolean). If true, use HTML table format.                   

% if varargin{k} = 'LaTeX', 
% then varargin{k+1} = logical (Boolean). If true, use LaTeX tabular format.        
                     
% if varargin{k} = 'Title', 
% then varargin{k+1} = character string = Title of table                                 

% if varargin{k} = 'ColName', 
% then varargin{k+1} = character string array = column headers
                                     
% if varargin{k} = 'File', 
% then varargin{k+1} = identifier of output file (default = screen)
                                       

%% INPUT PARSING
    
pa=inputParser; pa.CaseSensitive = false; % Initialize input parser

addRequired(pa,'a',                 @(x) isreal(x) && ismatrix(x));

addParameter(pa,'precision',[6],    @(x) all(isintegervalue(x)) && all(x)>0);
addParameter(pa,'integer',[],       @(x) all(isintegervalue(x)) && all(x)>0);
addParameter(pa,'normalized',false, @(x) islogical(x) && isscalar(x));
addParameter(pa,'startofline',' ',  @(x) ischar(x));
addParameter(pa,'separator',' ',    @(x) ischar(x));
addParameter(pa,'endofline',' ',    @(x) ischar(x));
addParameter(pa,'european',false,   @(x) islogical(x) && isscalar(x));
addParameter(pa,'csv',false,        @(x) islogical(x) && isscalar(x));
addParameter(pa,'html',false,       @(x) islogical(x) && isscalar(x));
addParameter(pa,'latex',false,      @(x) islogical(x) && isscalar(x));
addParameter(pa,'title',[],         @(x) ischar(x));
addParameter(pa,'colname',[],       @(x) iscellstr(x));
addParameter(pa,'file',[],          @(x) ischar(x));

parse(pa,a,varargin{:}); % Read the argument list
a=pa.Results.a;
p=pa.Results.precision;
intx=pa.Results.integer;
normalized=pa.Results.normalized;
sol=pa.Results.startofline;
sep=pa.Results.separator;
eol=pa.Results.endofline;

eur=pa.Results.european;
if eur
    euradix=' with European radix';
else
    euradix=[];
end

csv=pa.Results.csv;
if csv
    if eur
        sol=[]; sep=';'; eol=[];    
    else
        sol=[]; sep=','; eol=[];
    end
end

html=pa.Results.html;
if html
    sol='<tr align="right"><td>'; sep='</td><td>'; eol='</td></tr>';
end

latex=pa.Results.latex;
if latex
    sol=' '; sep=' & '; eol=' \\';
end

Title=pa.Results.title;
colname=pa.Results.colname;

file_name=pa.Results.file;
if isempty(file_name)
    fileID=1; % Reserved file handle for screen output
else
    fileID=fopen(file_name,'w');
end

%% CHECK INPUT CONSISTENCY

lp=length(p);

warning('off','backtrace')
for i=1:lp
    if (p(i)>15)
        warning(['The precision requested for printing ' ...
                 'column %d exceeds \n         '...
                 'the 15 significant digits that ' ...
                 'are warranted in double precision.'],i)
    end
    if (p(i)>17)
        warning(['The precision requested for printing ' ...
                 'column %d exceeds the maximum of \n         '...
                 '17 significant digits that ' ...
                 'can at best be achieved in double precision.'],i)
    end
end
warning('on','backtrace')

if (csv + html + latex) > 1
   error('CSV,HTML and LaTeX outputs can not be requested simultaneously.')
end

%% INITIALIZATIONS

[m,n]=size(a); 
l=sprintf('%d',ceil(log10(m)));

for j=1:n
    if ismember(j,intx)
        a(:,j)=round(a(:,j));
        k(j)=ceil(log10(max(abs(a(:,j)))))+1+length(sep); 
        if ~all(a(:,j)>=0)
            k(j)=k(j)+1;
        end
    else
        k(j)=0; % 
    end
end

if(lp>n)
    p(n+1:end)=[]; % Discard input in excess
end

if (lp<n)
    for i=lp+1:n
        p(i)=p(lp); % Complete if abbreviated input
    end
end

if normalized
    f='E'; % Normalized scientific notation forces decimal point after
           % first significant digit
    corr=-1;
else
    f='G'; % Scientific notation with floating position of decimal point
           % (default)
    corr=0;
end        

%% WIDTH OF DATA FIELDS

% For each column i,
% width(i) = Number of characters to print a floating-point number =
%  + Mantissa with p(i) significant digits: p(i)
%  + Sign of mantissa: 1  
%  + Decimal point:    1  
%  + Exponent flag:    1  
%  + Sign of exponent: 1  
%  + Exponent value:   3  
%  + character(s) separating fields > 0   

for i=1:n
    width(i)=p(i)+7+length(sep)+corr;
    s{i}=sprintf('%d',p(i)+corr); 
    w{i}=sprintf('%d',width(i)); 
end

%% BEGIN ENVIRONMENT

if (~isempty(Title) && ~csv && ~html && ~latex)
    fprintf(fileID,'\n%s\n',Title);
end

if csv
    fprintf('\nCSV formatting of table%s\n\n',euradix);
end

if html
    fprintf('\nHTML formatting of table%s\n\n',euradix);
    fprintf(fileID,'<table border="2" cellspacing="2" cellpadding="2">\n');
    fprintf(fileID,'<caption>%s</caption>\n',Title);
end

if latex
    fprintf('\nLaTeX formatting of table%s\n\n',euradix);
    fprintf(fileID,'\\begin{table}\n\\begin{tabular}{|');
    for j=1:n
        fprintf(fileID,'r|'); % Right alignment
    end
    fprintf(fileID,'}\n\\hline\n');
end

%% COLUMNS' HEADERS

if  isempty(colname) % No header
    
    KK=k; % Needed for formatting columns of integers
   
    for i=1:n % Padding with blanks not necessary
        D{i} =[]; DD{i}=[];
    end
    
else %  Headers activated by option 'ColName'
    
    if length(colname) ~= n
        for j=length(colname)+1:n
            colname{j}=' ';
        end
    end
    
    for j=1:n % Needed for formatting columns of integers
        KK(j)=max(k(j),length(colname{j}));
    end
    
    for j=1:n
        B=[];
        if width(j)>=length(colname{j})
            B=' ';
        end
        for i=1:(width(j)-1-length(colname{j}))
            B=[B ' '];
        end
        D{j}=B; % Blanks padding the field of the header of column j
                % to match the width of the data of column j
    end
    for j=1:n
        B=' ';
        if width(j)>=length(colname{j})
            B=[ ];
        end
        for i=1:(length(colname{j})-1-width(j))
            B=[B ' '];
        end
        DD{j}=B; % Blanks padding the field of the data of column j
                 % to match the length of the header of column j
    end
    
    % Print the line of column headers
    
    fprintf(fileID,'%s',sol);  
    for j=1:n-1
        if ismember(j,intx)
            B=[];
            for i=1:max(k(j),length(colname{j}))-length(colname{j})
                B=[B ' '];
            end
            fprintf(fileID,'%s%s%s',B,colname{j},sep);
        else
            fprintf(fileID,'%s%s%s',D{j},colname{j},sep);
        end
    end
    if ismember(n,intx)
        B=[];
        for i=1:max(k(j),length(colname{n}))-length(colname{n})
            B=[B ' '];
        end
        fprintf(fileID,'%s%s%s\n',B,colname{n},eol);
    else
        fprintf(fileID,'%s%s%s\n',D{n},colname{n},eol);
    end
    if latex
        fprintf(fileID,'\\hline\\hline\n');
    end
end

%% PRINT ROWS OF DATA

for i=1:m
    fprintf(fileID,'%s',sol);
    %% Standard columns 
    for j=1:n-1
        if ismember(j,intx) % Integer output
            w{j}=sprintf('%d',KK(j)); 
            fprintf(fileID,['%s' '%' w{j} '.0f' sep],DD{j},a(i,j));
        else % Output with fixed significant digit
            buf=sprintf(['%s' '%#' w{j} '.' s{j} f sep],DD{j},a(i,j));
            if eur
                buf=strrep(buf,'.',','); % Change radix
            end
            fprintf(fileID,buf);
        end
    end
    %% Last column features a specific termination
    if ismember(n,intx) % Integer output
        w{n}=sprintf('%d',KK(n)); 
        fprintf(fileID,['%s' '%' w{n} '.0f'],DD{n},a(i,n));
    else % Output with fixed significant digit 
        buf=sprintf(['%s' '%#' w{n} '.' s{n} f],DD{n},a(i,n));
        if eur
            buf=strrep(buf,'.',','); % Change radix
        end
        fprintf(fileID,buf);
    end
    fprintf(fileID,'%s \n',eol);
end

%% END OF ENVIRONMENT

if (~html && ~latex)
    fprintf(fileID,'\n');
end

if html
    fprintf(fileID,'</table>\n\n');
end

if latex
    fprintf(fileID,'\\hline\n\\end{tabular}\n\\begin{caption}\n');
    fprintf(fileID,'%s\n\\end{caption}\n\\end{table}\n\n',Title);
end

if (fileID ~=1)  
    fclose(fileID); 
end

end % End of function printtable

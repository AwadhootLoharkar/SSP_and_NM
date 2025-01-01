function savepdf(file_name,varargin)

%% PDF OUTPUT OF CURRENT FIGURE
%
%  WARNING: In order to work, a figure must have been 
%           opened by the calling script !

%% INPUT
%  file_name = file name with the structure <file_name_no_ext>.<ext>
%              NB: File names without any extension can be processed.
%            = file name that will be used to construct the output
%              file name as <file_name_no_ext>.tex
    
%% OUTPUT
%  <file_name_no_ext>.pdf file containing a figure to be used 
%  as follows in a LaTeX document: 

p=inputParser; p.CaseSensitive = false; % Initialize input parser
addRequired(p,'file_name', @(x) ischar(x));
parse(p,file_name); % Read the argument list
file_name=p.Results.file_name;

h = findall(0,'type','figure');
if isempty(h)
    error(['Can not produce pdf file ' ...
           'because no figure is open.'])
end

pdf_file=changefiletype(file_name,'pdf');

if isoctaveshell
    print(pdf_file,'-dpdf','-landscape','-fillpage');
end

if ismatlabshell
    exportgraphics(gca,pdf_file) % NB: Matlab introduced this function 
                                 % in 2020, to be preferred to
                                 %     print(pdf_file,'-dpdf',...)  
                                 % to enable a correct management of
                                 % computer fonts. Moreover, it
                                 % features a tight cropping (white
                                 % spaces around plot are minimized),
                                 % thereby making easier to include
                                 % the pdf file in a text document.
end

fprintf('Plot of window "%s" saved in %s\n',get(gcf,'name'),pdf_file);

end % End of function savepdf
function [new_file_name]=changefiletype(file_name,newftype)

%% CHANGING FILE TYPE IN THE STRING OF A FILE NAME
    
%% INPUT
% 
%  file_name = character string 
%            = file name structured as fname.ftype

%  newftype  = character string 
%            = targeted file type extension
    
%% OUTPUT
%   
%  new_file_name = character string 
%                = new file file name structured as fname.newftype
    
%  NB: The procedure works also if length(ftype)==0
%      In this case .newftype is appended to file_name=fname

%% INPUT PARSER
p=inputParser; p.CaseSensitive = false; % Initialize input parser
addRequired(p,'file_name', @(x) ischar(x));
addRequired(p,'newftype',  @(x) ischar(x));
parse(p,file_name,newftype); % Read the argument list
file_name=p.Results.file_name;
newftype=p.Results.newftype;

%% CORE JOB

dot=max(strfind(file_name,'.')); % Position of last dot in file_name
if isempty(dot)
    fname=file_name;  % If no file type given as input
else
    fname=substr(file_name,1,dot-1); % File name without filetype
                                     % extension
end
 
new_file_name=sprintf('%s.%s',fname,newftype); % New file name

end % End of function changefiletype

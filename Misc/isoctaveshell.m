function [OCTAVE]=isoctaveshell()
    OCTAVE=~isempty(strfind(evalc('ver'),'Octave')); % True if OCTAVE shell
end
% hello.m

mystartdefaults

char='Hello !';                % Character string variable
print_to_file= true;           % Logical (Boolean) Variable

fprintf('%s \n',char);         % Display  on screen

if (print_to_file)
  f1 = fopen ('out.txt','w');  % Open output file and
			       % attribute file identifier
  fprintf(f1,'%s \n',char);    % Print to output file
  fclose(f1);                  % Close output file
end

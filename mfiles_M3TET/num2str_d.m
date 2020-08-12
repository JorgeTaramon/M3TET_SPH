function string_number = num2str_d(number, digits)
%NUM2STR_D converts a number into a string with a given #of digits
% num2str_d(3, 5) will return the string '00003'
% This is mainly used in naming of files
%
% 25. January, 2000, Dani Schmid

string_number = num2str(number,['%',num2str(digits),'.',num2str(digits),'d']);

end
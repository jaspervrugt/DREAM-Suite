function strnumnew = addcomma(Num)
% Adapted by JAV to handle multiple numbers
strnumnew = cell(size(Num));
for zz = 1:numel(Num)
    num = Num(zz);
    
    if isnumeric(num)==1
        strnum=num2str(num);
        strnum2=regexprep(fliplr(strnum),' ','');
    else
        if ischar(num)==1
            if sum(ismember(num,','))==0
                strnum2=regexprep(fliplr(num),' ','');
            end
        end
    end
    if isempty(str2double(strnum2))==0
        strnum3=[];
        for i=1:size(strnum2,2)
            iv=logical(abs(double(logical(i-fix(i/3)*3))-1));
            if iv
                if i~=size(strnum2,2)
                    strv=[strnum2(1,i) ','];
                else
                    strv=strnum2(1,i);
                end
            else
                strv=strnum2(1,i);
            end
            strnum3 = [strnum3 strv];
        end
        strnumnew{zz} = fliplr(strnum3);
    else
        strnumnew{zz} = num;
    end
end
%end function

% Version that is much simpler and supposed to work
% % function numOut = addComma(numIn)
% %    jf=java.text.DecimalFormat; % comma for thousands, three decimal places
% %    numOut= char(jf.format(numIn)); % omit "char" if you want a string out
% % end
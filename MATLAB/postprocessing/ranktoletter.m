function c = letterRankToLetter(n)

d = ceil(n/26); 
c1 = repmat('A',1,d-1); nleft = n - (d-1)*26;
c2 = char('A'+nleft-1);
c = strcat(c1,c2);
% % switch d
% %     case 1
% %         c = char('A'+n-1);
% %     case 2
% %         
% % c = char('A'+n-1);
% % end

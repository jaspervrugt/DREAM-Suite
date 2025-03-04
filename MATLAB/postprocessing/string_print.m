function [ y ] = string_print ( x ); 
% Makes sure that plotting goes ok

n_char = 5; % for printing

for i = 1 : numel ( x ),

    % Extra space for minus or not?
    if x(i) < 0, 
        n_space = 1; 
    else
        n_space = 0;
    end
    % How many total places?
    n_tot = n_space + n_char;
    
    % Now get rid of e+ or e- operator
    a = sprintf('%6f',x(i));
    
    % Now do the job
    a = a(1:n_tot);
    
    % Return
    y(i) = {a};
    
end

% % Now loop over elements of x
% for i = 1 : numel ( x );
% 
%     % Check whether we need extra space or not
%     if x(i) < 0, 
%         n_space = 1; 
%     else
%         n_space = 0;
%     end
%     
%     % Now make string of number
%     p = num2str(x(i));
%     % How many elements does string have?
%     n_p = numel(p);
%     % How many elements total?
%     n_tot = n_char + n_space;
% 
%     % Now create extra zeros if not sufficient elements
%     if n_p < n_tot,
%         add = '0';
%         for j = 1 : n_tot - n_p - 1,
%             add = strcat(add,'0');
%         end
%         p = strcat(p,add);
%     else
%         p = p(1:n_tot);
%     end
%     % Now store in return argument
%     y(i) = {p};
% end
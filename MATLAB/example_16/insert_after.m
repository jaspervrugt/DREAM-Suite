function NewCharArray = insert_after( CharArray, Position, WhatToInsert)
% Inserts character into an existing character array

NewCharArray = char( strcat( cellstr(CharArray(:,1:Position)), ...
    cellstr(WhatToInsert), cellstr(CharArray(:, Position+1:end)) ) );
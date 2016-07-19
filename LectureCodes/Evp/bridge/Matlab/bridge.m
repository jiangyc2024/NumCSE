% Data for truss structure ``bridge''
pos = [ 0 0; 1 0;2 0;3 0;4 0;5 0;1 1;2 1;3 1;4 1;2.5 0.5];
con = [1 2;2 3;3 4;4 5;5 6;1 7;2 7;3 8;2 8;4 8;5 9;5 10;6 10;7 8;8 9;9 10;8 11 ...
       ; 9 11;3 11;4 11;4 9];
n = size(pos,1);
top = sparse(con(:,1),con(:,2),ones(size(con,1),1),n,n);
top = sign(top+top');

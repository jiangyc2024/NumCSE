function [val,row_ind,col_ptr] = CCS(A)
[row_ind col_ind val] = find(A);
col_ptr = find (diff ([0 ; col_ind ; size(A,2)+1]));

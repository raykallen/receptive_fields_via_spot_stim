function [neighb] = neighbormap( matrix )
%makes a neighbor map based off the size of the matrix

[n_rows, n_col] = size(matrix);
%neighb = zeros(n_rows, n_col);

for row = 1: n_rows
    for col = 1: n_col
        if row == 1 && col == 1 %corner
            neighb(1,1).map = [2, 1; 1, 2; 2, 2; row, col];
        elseif row == 1 && col == n_col %corner
                neighb(1, col).map = [1, (col-1); 2, (col -1); 2, col; row, col];
        elseif row == n_rows && col == 1 %corner
                neighb(row, 1).map = [(n_rows -1), 1; (n_rows -1), 2; n_rows, 2; row, col];
        elseif row == n_rows && col == n_col %corner
                neighb(row, col).map= [(row -1), (col -1); row, (col -1); (row - 1), col; row, col];
        elseif row == 1 %top row
                neighb(1, col).map = [1, (col -1); 1, (col +1); 2, (col-1); 2, col; 2, (col + 1); row, col];
        elseif col == 1 %far col
               neighb(row, 1).map = [(row-1), 1; (row+1), 1; (row-1), 2; row, 2; (row+1), 2; row, col];
        elseif row == n_rows %bottom row
               neighb(row, col).map = [row, (col-1); row, (col+1); (row-1), (col-1); (row-1), col; (row-1), (col+1); row, col];
        elseif col == n_col
               neighb(row,col).map = [(row-1), col; (row+1), col; (row-1), (col-1); row, (col-1); (row+1), (col-1); row, col];
        else
            neighb(row, col).map = [(row-1), (col-1); (row-1), col; (row-1), (col+1); row, (col-1); row, (col+1); (row+1), (col-1); (row+1), col; (row+1), (col+1); row, col];
        end
    end                            
end


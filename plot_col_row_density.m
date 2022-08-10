% asdf
load('program.mat')
A = A(1:v,1:v);
a = sort(full(sum(A != 0, 1)), 'descend');
b = sort(full(sum(A != 0, 2)), 'descend');



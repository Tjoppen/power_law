% Computes initial point by going from origin in the direction of c half as far as possible, then centering
function x = program_init(A, b, c)
  w = w_to_border(A, b, c, zeros(size(A,1),1));
  printf("Computing initial center\n");
  x = goto_center(A, b, w*c*0.5);
endfunction

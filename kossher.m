clear, clc
A = [1 1 -1 -1 0 0 0 0 0 0 0;
     0 1 0 0 0 0 0 0 0 0 0;
     0 0 1 0 -1 0 -1 0 0 0 0;
     0 0 0 1 0 -1 0 -1 0 0 0;
     0 0 0 0 1 0 0 0 0 0 0;
     0 0 0 0 0 1 0 0 0 0 0;
     0 0 0 0 0 0 1 0 0 0 0;
     0 0 0 0 0 0 0 1 0 0 0;
     0 0 0 0 0 0 0 0 1 0 0;
     0 0 0 0 0 0 -1 -1 -1 1 0;
     0 0 0 0 0 0 0 0 0 0 1;];


c = [103; 52; 13+5*9; 33+57; 17; 33; 13+5*9-17; 57; 41; 117; 117;];


linAns = linsolve(A,c);

ASE = sum((linAns-c).^2)
RSE = sum(((linAns-c)./linAns).^2)


% Define the process matrix A and measurement vector c
A = [1  1 -1 -1  0  0  0  0  0  0  0;
     0  1  0  0  0  0  0  0  0  0  0;
     0  0  1  0 -1  0 -1  0  0  0  0;
     0  0  0  1  0 -1  0 -1  0  0  0;
     0  0  0  0  1  0  0  0  0  0  0;
     0  0  0  0  0  1  0  0  0  0  0;
     0  0  0  0  0  0  1  0  0  0  0;
     0  0  0  0  0  0  0  1  0  0  0;
     0  0  0  0  0  0  0  0  1  0  0;
     0  0  0  0  0  0 -1 -1 -1  1  0;
     0  0  0  0  0  0  0  0  0  0  1;];

c = [103; 52; 58; 90; 17; 33; 41; 57; 41; 117; 117;];

% Define the objective function (least squares)
W = eye(length(c)); % Assuming equal weights for simplicity
f_R = fmincon(@(f_R) (f_R - c)' * W * (f_R - c), c, [], [], A, zeros(size(A, 1), 1));

% Display the reconciled values
disp('Reconciled values:');
disp(f_R);

% EXAMPLE_MICROSTATE_SURROGATES
% Minimal example for generating S1-S7 surrogate microstate sequences.

seq = [1 1 1 2 2 3 3 3 3 1 1 4 4 2 2 2 3 3 1 1 1 1 4 4 4 2 2];
states = 1:4;
seed = 1;
Lm = 3;

S1 = generate_microstate_surrogate(seq, 'S1', 'Lm', Lm);
S2 = generate_microstate_surrogate(seq, 'S2', 'states', states, 'seed', seed);
S3 = generate_microstate_surrogate(seq, 'S3', 'states', states, 'seed', seed);
S4 = generate_microstate_surrogate(seq, 'S4', 'seed', seed, 'avoidAdjacent', false);
S5 = generate_microstate_surrogate(seq, 'S5', 'Lm', Lm, 'seed', seed);
S6 = generate_microstate_surrogate(seq, 'S6', 'states', states, 'seed', seed);
S7 = generate_microstate_surrogate(seq, 'S7', 'states', states, 'seed', seed);

disp('Original:'); disp(seq);
disp('S1:'); disp(S1);
disp('S2:'); disp(S2);
disp('S3:'); disp(S3);
disp('S4:'); disp(S4);
disp('S5:'); disp(S5);
disp('S6:'); disp(S6);
disp('S7:'); disp(S7);

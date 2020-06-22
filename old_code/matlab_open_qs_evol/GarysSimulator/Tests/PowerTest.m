propagator.precalculation(params.Trotter);
tic
u=propagator.unitary();
whos u
spy(u{1})
toc
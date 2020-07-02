%% Total variation regularisation in tomography applications
% This is a little program that allows to reconstruct a 3d image from a set
% of 2d projections. The program therefore realises the minimisation of the
% mathematical model
%
% $$ u \in \arg\min_{u \geq 0} \left\{ \frac{\lambda}{2}\| Ru - f \|_2^2 +
% \| \| \nabla u \|_2 \|_1 \right\} $$
%
% via a preconditioned primal dual hybrid gradient method (PDHGM) as
% described in [1]. In order to use the algorithm described in [1], we have
% to formulate our mathematical problem as
%
% $$ \min F(Kx) + G(x) $$
%
% with $K = (R, \nabla)^T$, $F(z_1, z_2) = \lambda \| z_1 - f \|_2^2 / 2 +
% \| \| z_2 \|_2 \|_1$, and $G(x) = \chi_P(x)$. Here $\chi_P$ denotes the
% characteristic function over the set $P := \{ x \ | \ x \geq 0 \}$. The
% PDHGM does not compute a solution of the primal problem, but of the
% primal-dual problem
%
% $$ \min \max \langle Kx, y \rangle + G(x) - F^*(y) ;$$
%
% here $F^*$ is the Fenchel-dual of $F$ which is in our case given as
% $F^*(y_1, y_2) = \| y_1 \|_2^2 / (2\lambda) + \langle y_1, f\rangle +
% \chi_U(y_2)$, where $\chi_U$ denotes the characteristic function over the
% set $U := \{ y \ | \ \| \| y \|_2 \|_\infty \leq 1\}$. Given the matrix
% $K$ and the resolvent operators of $F^*$ and $G$, we compute the iterates
% as described in [1, Equation (4)].
%
% Author: Martin Benning - mb941@cam.ac.uk
%
% Date: 23.12.15
%
%% Load projection images and create radon operator
% We are going to load exemplary projection images that were made from the
% MATLAB? dataset *wmri*. We further reduce the number of sampled angles to
% 45; we only use every fourth projection of the original dataset that
% included 180 projections in total (for the angular range 0 to 179).
load('projectionimages.mat')
angles = (0:4:179)';
imsize = [128, 128, size(proj, 2)];
projsmall = proj(:, :, 1:4:180);
%%
% Now we compute a matrix that maps the 3d image onto the sinogram of each
% z-component. It is recommended to open the parallel pool before executing
% the following line of code.
tic; R = radonmtx(imsize(1), angles, size(proj, 1)); toc;
%%
% Subsequently, we initalise the operator that maps the 3d image onto the
% projection images, and store the projection images as our data *f* in
% terms of a column vector.
Rop = matleftmult(R, [size(R, 2) size(proj, 2)]);
f = permute(projsmall, [1 3 2]);
f = f(:); %sinogram data for each z-slice
scal = 4096; %rescaling of the data; acts as a regularisation parameter
%% Initialise 3d gradient operator
% We initialise a 3d forward finite-difference approximation of the
% gradient operator with the following command.
Grad = fcthdlop(imsize, [imsize 3], @fwgradient3d, @bwgradient3d);
%% Initialise operator K
% In order to initialise the operator $K$, we only have to concatenate the
% operators *Rop* and *Grad*.
K = [Rop; Grad];
%% Initialise operators for F* and G
% Given the operator $K$, we need to initalise objects of the functionals
% $F^*$ and $G$ that allow the computation of the
% resolvent/proximity-operations.
Fstar = dualvpproj(imsize, size(projsmall));
Fstar.setproxdata(f*scal)
G = nonnegproj(imsize);
%% Initialise PDHGM
% Now we have all the ingredients to initalise an instance of the PDHGM. We
% set the number of iterations arbitrarily to 300.
solver = pdhgm(K, Fstar, G);
solver.setmaxiter(300)
%% Compute preconditioners
% Before we execute the solver, we need to choose parameters $T$ and
% $\Sigma$, as they are required for the execution of [1, Equation (4)]. We
% have precomputed those in analogy to the example in [1, Lemma 2], for
% $\alpha = 1$, and stored the parameters in 'sigmatau.mat'. We load these
% parameters and pass them to the objects representing $F^*$ and $G$.
load('sigmatau.mat')
Fstar.setproxparam(sigma1, sigma2);
G.setproxparam(tau);
%% Run PDHGM
% Now we run 300 iterations of the PDHGM. We could change the value of the
% regulrisation parameter $\lambda$; however, in this setup the re-scaling
% of the initial data $f$ acts also as a regularisation. As we found the
% iteration to be more stable with $\lambda = 1$, we stick to the
% re-scaling of $f$.
lambda = 1;
Fstar.setregularisationparameter(lambda)
solver.solve
%% Visualise results
% To conlcude, we visuale the 13-th slice of the variable $x$ after 300
% iterations of PDHGM.
x = reshape(solver.getvariables.x, imsize);
imagesc(x(:, :, 13))
axis image
colorbar
colormap(gray(512))
drawnow
%% References
% [1] Pock, Thomas, and Antonin Chambolle. "Diagonal preconditioning for
% first order primal-dual algorithms in convex optimization." Computer
% Vision (ICCV), 2011 IEEE International Conference on. IEEE, 2011.
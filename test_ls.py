#Required packages
#%%
import os
import sys
sys.path.append('./src')
from supportfunctions import *
sys.stdout.flush()
import petsc4py
petsc4py.init(sys.argv)
from petsc4py import PETSc
from scipy.sparse import spdiags
from scipy.sparse import coo_matrix
from scipy.sparse import csr_matrix
from datetime import datetime
# Damage function choices
damageSpec = 'Weighted'  # Choose among "High"(Weitzman), 'Low'(Nordhaus) and 'Weighted' (arithmeticAverage of the two)

if damageSpec == 'High':
    weight = 0.0
elif damageSpec == 'Low':
    weight = 1.0
else:
    weight = 0.5

ξp =  1 / 4000  # Ambiguity Averse Paramter
# We stored solutions for ξp =  1 / 4000 to which we referred as "Ambiguity Averse" and ξp = 1000 as “Ambiguity Neutral” in the paper
# Sensible choices are from 0.0002 to 4000, while for parameters input over 0.01 the final results won't alter as much

if ξp < 1:
    aversespec = "Averse"
else:
    aversespec = 'Neutral'

smart_guess = False
model = 'petsc'
current_time = datetime.now()
filename = 'BBH_' + model + '_' + damageSpec + '_' + aversespec + '_' + "{}-{}-{}".format(current_time.day, current_time.hour, current_time.minute)

McD = np.loadtxt('./data/TCRE_MacDougallEtAl2017_update.txt')
par_lambda_McD = McD / 1000

β𝘧 = np.mean(par_lambda_McD)  # Climate sensitivity parameter, MacDougall (2017)
σᵦ = np.var(par_lambda_McD, ddof = 1)  # varaiance of climate sensitivity parameters
λ = 1.0 / σᵦ

quadrature = 'legendre'
n = 30
a = β𝘧 - 5 * np.sqrt(σᵦ)
b = β𝘧 + 5 * np.sqrt(σᵦ)

(xs,ws) = quad_points_legendre(n)
xs = (b-a) * 0.5  * xs + (a + b) * 0.5
s = np.prod((b-a) * 0.5)

start_time = time.time()

# Parameters as defined in the paper
δ = 0.01
κ = 0.032
σ𝘨 = 0.02
σ𝘬 = 0.0161
σ𝘳 = 0.0339
α = 0.115000000000000
ϕ0 = 0.0600
ϕ1 = 16.666666666666668
μk = -0.034977443912449
ψ0 = 0.112733407891680
ψ1 = 0.142857142857143

# parameters for damage function settings
power = 2
γ1 = 0.00017675
γ2 = 2. * 0.0022
γ2_plus = 2. * 0.0197
γ̄2_plus = weight * 0 + (1 - weight) * γ2_plus

σ1 = 0
σ2 = 0
ρ12 = 0
F̄ = 2
crit = 2
F0 = 1

xi_d = -1 * (1 - κ)

# See Remark 2.1.1 regarding the choice of ε and η
# False Trasient Time step
ε = 0.3
# Cobweb learning rate
η = 0.05


# Specifying Tolerance level
tol = 1e-8

# Grids Specification

# Coarse Grids
# R_min = 0
# R_max = 9
# F_min = 0
# F_max = 4000
# K_min = 0
# K_max = 18
# nR = 4
# nF = 4
# nK = 4
# R = np.linspace(R_min, R_max, nR)
# F = np.linspace(F_min, F_max, nF)
# K = np.linspace(K_min, K_max, nK)

# hR = R[1] - R[0]
# hF = F[1] - F[0]
# hK = K[1] - K[0]

# Dense Grids

R_min = 0
R_max = 9
F_min = 0
F_max = 4000
K_min = 0
K_max = 18

hR = 0.05
hF = 25
hK = 0.15

R = np.arange(R_min, R_max + hR, hR)
nR = len(R)
F = np.arange(F_min, F_max + hF, hF)
nF = len(F)
K = np.arange(K_min, K_max + hK, hK)
nK = len(K)

# Discretization of the state space for numerical PDE solution.
# See Remark 2.1.2
(R_mat, F_mat, K_mat) = np.meshgrid(R,F,K, indexing = 'ij')
stateSpace = np.hstack([R_mat.reshape(-1,1,order = 'F'),F_mat.reshape(-1,1,order = 'F'),K_mat.reshape(-1,1,order = 'F')])

# Inputs for function quad_int
# Integrating across parameter distribution
quadrature = 'legendre'
n = 30
a = β𝘧 - 5 * np.sqrt(σᵦ)
b = β𝘧 + 5 * np.sqrt(σᵦ)

v0 = κ * R_mat + (1-κ) * K_mat - β𝘧 * F_mat

FC_Err = 1
episode = 0

if smart_guess:
    v0 = v0_guess
    q = q_guess
    e_star = e_guess
    episode = 1
    ε = 0.2

#%%
while FC_Err > tol and episode < 2:
    vold = v0.copy()
    # Applying finite difference scheme to the value function
    v0_dr = finiteDiff(v0,0,1,hR)
    v0_df = finiteDiff(v0,1,1,hF)
    v0_dk = finiteDiff(v0,2,1,hK)

    v0_drr = finiteDiff(v0,0,2,hR)
    v0_drr[v0_dr < 1e-16] = 0
    v0_dr[v0_dr < 1e-16] = 1e-16
    v0_dff = finiteDiff(v0,1,2,hF)
    v0_dkk = finiteDiff(v0,2,2,hK)
    if episode > 2000:
        ε = 0.1
    elif episode > 1000:
        ε = 0.2
    else:
        pass

    if episode == 0:
        # First time into the loop
        B1 = v0_dr - xi_d * (γ1 + γ2 * F_mat * β𝘧 + γ2_plus * (F_mat * β𝘧 - F̄) ** (power - 1) * (F_mat >= (crit / β𝘧))) * β𝘧 * np.exp(R_mat) - v0_df * np.exp(R_mat)
        C1 = - δ * κ
        e = -C1 / B1
        e_hat = e
        Acoeff = np.ones(R_mat.shape)
        Bcoeff = ((δ * (1 - κ) * ϕ1 + ϕ0 * ϕ1 * v0_dk) * δ * (1 - κ) / (v0_dr * ψ0 * 0.5) * np.exp(0.5 * (R_mat - K_mat))) / (δ * (1 - κ) * ϕ1)
        Ccoeff = -α  - 1 / ϕ1
        j = ((-Bcoeff + np.sqrt(Bcoeff ** 2 - 4 * Acoeff * Ccoeff)) / (2 * Acoeff)) ** 2
        i = α - j - (δ * (1 - κ)) / (v0_dr * ψ0 * 0.5) * j ** 0.5 * np.exp(0.5 * (R_mat - K_mat))
        q = δ * (1 - κ) / (α - i - j)
    else:
        e_hat = e_star

        # Step 4 (a) : Cobeweb scheme to update controls i and j; q is an intermediary variable that determines i and j
        Converged = 0
        nums = 0
        while Converged == 0:
            i_star = (ϕ0 * ϕ1 * v0_dk / q - 1) / ϕ1
            j_star = (q * np.exp(ψ1 * (R_mat - K_mat)) / (v0_dr * ψ0 * ψ1)) ** (1 / (ψ1 - 1))
            if α > np.max(i_star + j_star):
                q_star = η * δ * (1 - κ) / (α - i_star - j_star) + (1 - η) * q
            else:
                q_star = 2 * q
            if np.max(abs(q - q_star) / η) <= 1e-5:
                Converged = 1
                q = q_star
                i = i_star
                j = j_star
            else:
                q = q_star
                i = i_star
                j = j_star

            nums += 1
        if episode % 100 == 0:
            print('Cobweb Passed, iterations: {}, i error: {:10f}, j error: {:10f}'.format(nums, np.max(i - i_star), np.max(j - j_star)))

    i[i <= -1/ϕ1] = - 1/ϕ1 + 1e-8

    a1 = np.zeros(R_mat.shape)
    b1 = xi_d * e_hat * np.exp(R_mat) * γ1
    c1 = 2 * xi_d * e_hat * np.exp(R_mat) * F_mat * γ2
    λ̃1 = λ + c1 / ξp
    β̃1 = β𝘧 - c1 * β𝘧 / (ξp * λ̃1) -  b1 /  (ξp * λ̃1)
    I1 = a1 - 0.5 * np.log(λ) * ξp + 0.5 * np.log(λ̃1) * ξp + 0.5 * λ * β𝘧 ** 2 * ξp - 0.5 * λ̃1 * (β̃1) ** 2 * ξp
    R1 = 1 / ξp * (I1 - (a1 + b1 * β̃1 + c1 / 2 * β̃1 ** 2 + c1 / 2 / λ̃1))
    J1_without_e = xi_d * (γ1 * β̃1 + γ2 * F_mat * (β̃1 ** 2 + 1 / λ̃1)) * np.exp(R_mat)

    π̃1 = weight * np.exp(-1 / ξp * I1)

    # Step (2), solve minimization problem in HJB and calculate drift distortion
    # See remark 2.1.3 for more details
    start_time2 = time.time()
    if episode == 0 or (smart_guess and episode == 1):
        #@nb.jit(nopython = True, parallel = True)
        def scale_2_fnc(x, ndist, e_hat):
            return np.exp(-1 / ξp * x * e_hat) * ndist

        #@nb.jit(nopython = True, parallel = True)
        def q2_tilde_fnc(x, e_hat, scale_2):
            return np.exp(-1 / ξp * x * e_hat) / scale_2

        #@nb.jit(nopython = True, parallel = True)
        def J2_without_e_fnc(x, ndist, e_hat, scale_2):
            return x * q2_tilde_fnc(x, e_hat, scale_2) * ndist

        (xs,ws) = quad_points_legendre(n)
        xs = (b-a) * 0.5  * xs + (a + b) * 0.5
        s = np.prod((b-a) * 0.5)

        normdists = np.zeros(n)
        distort_terms = np.zeros((n, nR, nF, nK))
        for i_iter in range(n):
            normdists[i_iter] = normpdf(xs[i_iter],β𝘧,np.sqrt(σᵦ))
            distort_terms[i_iter] = xi_d * (γ1 * xs[i_iter] + γ2 * xs[i_iter] ** 2 * F_mat + γ2_plus * xs[i_iter] * (xs[i_iter] * F_mat - F̄) ** (power - 1) * ((xs[i_iter] * F_mat - F̄) >= 0)) * np.exp(R_mat)

        dVec = [hR, hF, hK]
        increVec = [1, nR, nR * nF]
        I_LB_R = (stateSpace[:,0] == R_min)
        I_UB_R = (stateSpace[:,0] == R_max)
        I_LB_F = (stateSpace[:,1] == F_min)
        I_UB_F = (stateSpace[:,1] == F_max)
        I_LB_K = (stateSpace[:,2] == K_min)
        I_UB_K = (stateSpace[:,2] == K_max)

    scale_2 = np.zeros(F_mat.shape)
    for i_iter in range(n):
        scale_2 += ws[i_iter] * scale_2_fnc(distort_terms[i_iter], normdists[i_iter], e_hat)
    scale_2 = s * scale_2

    I2 = -1 * ξp * np.log(scale_2)

    J2_without_e = np.zeros(F_mat.shape)
    for i_iter in range(n):
        J2_without_e += ws[i_iter] * J2_without_e_fnc(distort_terms[i_iter], normdists[i_iter], e_hat, scale_2)
    J2_without_e = s * J2_without_e
    J2_with_e = J2_without_e * e_hat
    end_time2 = time.time()


    R2 = (I2 - J2_with_e) / ξp
    π̃2 = (1 - weight) * np.exp(-1 / ξp * I2)
    π̃1_norm = π̃1 / (π̃1 + π̃2)
    π̃2_norm = 1 - π̃1_norm

    # step 4 (b) updating e based on first order conditions
    expec_e_sum = (π̃1_norm * J1_without_e + π̃2_norm * J2_without_e)

    B1 = v0_dr - v0_df * np.exp(R_mat) - expec_e_sum
    C1 = -δ * κ
    e = -C1 / B1
    e_star = e

    J1 = J1_without_e * e_star
    J2 = J2_without_e * e_star

    # Step (3) calculating implied entropies
    I_term = -1 * ξp * np.log(π̃1 + π̃2)

    R1 = (I1 - J1) / ξp
    R2 = (I2 - J2) / ξp

    # Step (5) solving for adjusted drift
    drift_distort = (π̃1_norm * J1 + π̃2_norm * J2)

    if weight == 0 or weight == 1:
        RE = π̃1_norm * R1 + π̃2_norm * R2
    else:
        RE = π̃1_norm * R1 + π̃2_norm * R2 + π̃1_norm * np.log(
            π̃1_norm / weight) + π̃2_norm * np.log(π̃2_norm / (1 - weight))

    RE_total = ξp * RE

    # Step (6) and (7) Formulating HJB False Transient parameters
    # See remark 2.1.4 for more details
    A = -δ * np.ones(R_mat.shape)
    B_r = -e_star + ψ0 * (j ** ψ1) * np.exp(ψ1 * (K_mat - R_mat)) - 0.5 * (σ𝘳 ** 2)
    B_f = e_star * np.exp(R_mat)
    B_k = μk + ϕ0 * np.log(1 + i * ϕ1) - 0.5 * (σ𝘬 ** 2)
    C_rr = 0.5 * σ𝘳 ** 2 * np.ones(R_mat.shape)
    C_ff = np.zeros(R_mat.shape)

    C_kk = 0.5 * σ𝘬 ** 2 * np.ones(R_mat.shape)
    D = δ * κ * np.log(e_star) + δ * κ * R_mat + δ * (1 - κ) * (np.log(α - i - j) + K_mat) + drift_distort + RE_total # + I_term



    out_cpp = PDESolver(stateSpace, A, B_r, B_f, B_k, C_rr, C_ff, C_kk, D, v0, ε, solverType = 'False Transient')
    out_comp_cpp = out_cpp[2].reshape(v0.shape,order = "F")



    # Transforming the 3-d coefficient matrix to 1-dimensional
    A = A.reshape(-1,1,order = 'F')
    B = np.hstack([B_r.reshape(-1,1,order = 'F'),B_f.reshape(-1,1,order = 'F'),B_k.reshape(-1,1,order = 'F')])
    C = np.hstack([C_rr.reshape(-1,1,order = 'F'), C_ff.reshape(-1,1,order = 'F'), C_kk.reshape(-1,1,order = 'F')])
    D = D.reshape(-1,1,order = 'F')
    v0 = v0.reshape(-1,1,order = 'F')


    # prepare for ksp
    B_plus = np.maximum(B, np.zeros(B.shape))
    B_minus = np.minimum(B, np.zeros(B.shape))

    diag_0 = (A[:,0] - 1 / ε
              + (I_LB_R * B[:,0] / -dVec[0] + I_UB_R * B[:,0] / dVec[0] - (1 - I_LB_R - I_UB_R) * (B_plus[:,0] - B_minus[:,0]) / dVec[0] + (I_LB_R * C[:,0] + I_UB_R * C[:,0] - 2 * (1 - I_LB_R - I_UB_R) * C[:,0]) / dVec[0] ** 2)
              + (I_LB_F * B[:,1] / -dVec[1] + I_UB_F * B[:,1] / dVec[1] - (1 - I_LB_F - I_UB_F) * (B_plus[:,1] - B_minus[:,1]) / dVec[1] + (I_LB_F * C[:,1] + I_UB_F * C[:,1] - 2 * (1 - I_LB_F - I_UB_F) * C[:,1]) / dVec[1] ** 2)
              + (I_LB_K * B[:,2] / -dVec[2] + I_UB_K * B[:,2] / dVec[2] - (1 - I_LB_K - I_UB_K) * (B_plus[:,2] - B_minus[:,2]) / dVec[2] + (I_LB_K * C[:,2] + I_UB_K * C[:,2] - 2 * (1 - I_LB_K - I_UB_K) * C[:,2]) / dVec[2] ** 2))
    diag_R = (I_LB_R * B[:,0] / dVec[0] + (1 - I_LB_R - I_UB_R) * B_plus[:,0] / dVec[0] - 2 * I_LB_R * C[:,0] / dVec[0] ** 2 + (1 - I_LB_R - I_UB_R) * C[:,0] / dVec[0] ** 2)
    diag_Rm = (I_UB_R * B[:,0] / -dVec[0] - (1 - I_LB_R - I_UB_R) * B_minus[:,0] / dVec[0] - 2 * I_UB_R * C[:,0] / dVec[0] ** 2 + (1 - I_LB_R - I_UB_R) * C[:,0] / dVec[0] ** 2)
    diag_F = (I_LB_F * B[:,1] / dVec[1] + (1 - I_LB_F - I_UB_F) * B_plus[:,1] / dVec[1] - 2 * I_LB_F * C[:,1] / dVec[1] ** 2 + (1 - I_LB_F - I_UB_F) * C[:,1] / dVec[1] ** 2)
    diag_Fm = (I_UB_F * B[:,1] / -dVec[1] - (1 - I_LB_F - I_UB_F) * B_minus[:,1] / dVec[1] - 2 * I_UB_F * C[:,1] / dVec[1] ** 2 + (1 - I_LB_F - I_UB_F) * C[:,1] / dVec[1] ** 2)
    diag_K = (I_LB_K * B[:,2] / dVec[2] + (1 - I_LB_K - I_UB_K) * B_plus[:,2] / dVec[2] - 2 * I_LB_K * C[:,2] / dVec[2] ** 2 + (1 - I_LB_K - I_UB_K) * C[:,2] / dVec[2] ** 2)
    diag_Km = (I_UB_K * B[:,2] / -dVec[2] - (1 - I_LB_K - I_UB_K) * B_minus[:,2] / dVec[2] - 2 * I_UB_K * C[:,2] / dVec[2] ** 2 + (1 - I_LB_K - I_UB_K) * C[:,2] / dVec[2] ** 2)
    diag_RR = I_LB_R * C[:,0] / dVec[0] ** 2
    diag_RRm = I_UB_R * C[:,0] / dVec[0] ** 2
    diag_FF = I_LB_F * C[:,1] / dVec[1] ** 2
    diag_FFm = I_UB_F * C[:,1] / dVec[1] ** 2
    diag_KK = I_LB_K * C[:,2] / dVec[2] ** 2
    diag_KKm = I_UB_K * C[:,2] / dVec[2] ** 2

    #data = np.vstack([diag_0, diag_R, diag_Rm, diag_RR, diag_RRm, diag_F, diag_Fm, diag_FF, diag_FFm, diag_K, diag_Km, diag_KK, diag_KKm])
    data = [diag_0, diag_R, diag_Rm, diag_RR, diag_RRm, diag_F, diag_Fm, diag_FF, diag_FFm, diag_K, diag_Km, diag_KK, diag_KKm]
    diags = np.array([0,-increVec[0],increVec[0],-2*increVec[0],2*increVec[0],
                     -increVec[1],increVec[1],-2*increVec[1],2*increVec[1],
                     -increVec[2],increVec[2],-2*increVec[2],2*increVec[2]])

    # the transpose of matrix A_sp is the desired A matrix
    A_sp = spdiags(data, diags, len(diag_0), len(diag_0))
    b = -v0 / ε - D

    csr_mat = csr_matrix(A_sp.T)
    petsc_mat = PETSc.Mat().createAIJ(size=csr_mat.shape, csr=(csr_mat.indptr, csr_mat.indices, csr_mat.data))
    petsc_rhs = PETSc.Vec().createWithArray(b)
    x = petsc_mat.createVecRight()
    x.set(0)

    pickle.dump(csr_mat, open( "csr_mat", "wb"))
    np.save("b", b)
    # create linear solver
    start = time.time()
    ksp = PETSc.KSP()
    ksp.create(PETSc.COMM_WORLD)
    ksp.setType('lsqr')
    ksp.getPC().setType('jacobi')
    ksp.setOperators(petsc_mat)
    ksp.setFromOptions()
    ksp.setTolerances(rtol=1e-10, atol=1e-20)
    # petsc_time1 = time.time()
    ksp.solve(petsc_rhs, x)
    print(ksp.getConvergedReason())
    V = np.array(ksp.getSolution())
    out_comp = V.reshape(R_mat.shape,order = "F")
    A = A.reshape(R_mat.shape,order = "F")
    D = D.reshape(R_mat.shape,order = "F")
    v0 = v0.reshape(R_mat.shape,order = "F")
    # Calculating PDE error and False Transient error
    PDE_rhs = A * v0 + B_r * v0_dr + B_f * v0_df + B_k * v0_dk + C_rr * v0_drr + C_kk * v0_dkk + C_ff * v0_dff + D
    PDE_Err = np.max(abs(PDE_rhs))
    FC_Err = np.max(abs((out_comp - v0)))

    v_cpp = np.array(out_cpp[2])
    np.save("v_cpp", v_cpp)

    res_1 = np.linalg.norm(csr_mat@v_cpp - b)

    print(res_1)

    if episode % 1 == 0:
        print("Episode {:d}: PDE Error: {:.10f}; False Transient Error: {:.10f};" .format(episode, PDE_Err, FC_Err))
    print("Time for this iteration: {}".format(time.time() - start))
    episode += 1


    # step 9: keep iterating until convergence
    v0 = out_comp

print("Episode {:d}: PDE Error: {:.10f}; False Transient Error: {:.10f}" .format(episode, PDE_Err, FC_Err))
print("--- %s seconds ---" % (time.time() - start_time))

import pickle
# filename = filename
my_shelf = {}
for key in dir():
    if isinstance(globals()[key], (int,float, np.float, str, bool, np.ndarray,list)):
        try:
            my_shelf[key] = globals()[key]
        except TypeError:
            #
            # __builtins__, my_shelf, and imported modules can not be shelved.
            #
            print('ERROR shelving: {0}'.format(key))
    else:
        pass


file = open(filename, 'wb')
pickle.dump(my_shelf, file)
file.close()

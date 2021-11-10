import numpy as np
import scipy.linalg as slg
from numpy import tensordot as td
import qutip as qt
import matplotlib.pyplot as plt
import scipy.io as sio
import time

t_i = time.time()
Dmc = 5  # 腔模的维度
N, m = 105, 35
N_p = int(N / m)
Dms = N_p + 1
Ds2 = Dms ** 2
I_e = np.eye(Dms)
I_c = np.eye(Dmc)
jz = np.array(qt.jmat(N_p / 2, 'z'))
jp = np.array(qt.jmat(N_p / 2, '+'))
jm = np.array(qt.jmat(N_p / 2, '-'))
a = np.array(qt.destroy(Dmc))
II_e = td(I_e, I_e, 0)
II_c = td(I_c, I_c, 0)


def Generate_Exp_Lindblads(tau, N_prime, drive, Domega, G, gammah, kappa):
    """
    :param tau: 时间演化步长
    :param drive: 激光驱动强度
    :param Domega: 原子与激光的失谐，假设激光频率与腔频共振
    :param G: 原子与腔的耦合强度
    :param gamma: 原子耗散率
    :param kappa: 腔场耗散率
    :return: Lindblad演化算符列表
    """
    jl = [td(jp, I_e, 0), td(jm, I_e, 0), td(jz, I_e, 0)]
    jr = [td(I_e, jp.T, 0), td(I_e, jm.T, 0), td(I_e, jz.T, 0)]
    cl = [td(a.T, I_c, 0), td(a, I_c, 0)]
    cr = [td(I_c, a, 0), td(I_c, a.T, 0)]
    m = len(Domega)
    Exp_Lindblad = []
    for k in range(m):
        Lindblad = -1j * (Domega[k] * td(td(jz, I_e, 0), II_c, 0) + G[k] * (td(jl[0], cl[1], 0) + td(jl[1], cl[0], 0)) +
                          drive / m * td(II_e, cl[0] + cl[1], 0) -
                          (Domega[k] * td(td(I_e, jz, 0), II_c, 0) + G[k] * (td(jr[0], cr[1], 0) + td(jr[1], cr[0], 0)) +
                           drive / m * td(II_e, cr[0] + cr[1], 0))) + gammah / (2 * N_prime) * (
                           2 * td(td(jm, jm, 0), II_c, 0) - td(td(jp.dot(jm), I_e, 0), II_c, 0) -
                           td(td(I_e, jp.dot(jm), 0), II_c, 0)) + kappa / 2 * (
                           2 * td(II_e, td(a, a, 0), 0) - td(II_e, td(a.T.dot(a), I_c, 0), 0) -
                           td(II_e, td(I_c, a.T.dot(a), 0), 0)) / m
        Lindblad = Lindblad.transpose(0, 2, 4, 6, 1, 3, 5, 7).reshape(Ds2 * 25, Ds2 * 25)
        if k == 0:
            Exp_Lindblad.append(slg.expm(tau * Lindblad).T)
        else:
            Exp_Lindblad.append(slg.expm(tau / 2 * Lindblad).T)
        # for k in range(m):
        #     Lindblad = -1j * (Domega[k] * td(td(sp.dot(sm), I_e, 0), II_c, 0) + G[k] * (td(sl[0], cl[1], 0) + td(sl[1], cl[0], 0)) +
        #                       drive / (m * N_prime) * td(II_e, cl[0] + cl[1], 0) -
        #                       (Domega[k] * td(td(I_e, sp.dot(sm), 0), II_c, 0) + G[k] * (td(sr[0], cr[1], 0) + td(sr[1], cr[0], 0)) +
        #                        drive / (m * N_prime) * td(II_e, cr[0] + cr[1], 0))) + gammah / 2 * (
        #                        2 * td(td(sm, sm, 0), II_c, 0) - td(td(sp.dot(sm), I_e, 0), II_c, 0) -
        #                        td(td(I_e, sp.dot(sm), 0), II_c, 0)) + gammap * (td(td(sz, sz, 0), II_c, 0) - td(II_e, II_c, 0)) + kappa / 2 * (
        #                        2 * td(II_e, td(a, a, 0), 0) - td(II_e, td(a.T.dot(a), I_c, 0), 0) -
        #                        td(II_e, td(I_c, a.T.dot(a), 0), 0)) / (m * N_prime)
        # Lindblad = Lindblad.transpose(0, 2, 4, 6, 1, 3, 5, 7).reshape(100, 100)
        # if k == 0:
        #     Exp_Lindblad_1 = slg.expm(tau * Lindblad).T
        # Exp_Lindblad.append(slg.expm(tau / 2 * Lindblad).T)
    return Exp_Lindblad


def update_coff(L, R, C):
    return np.squeeze(td(L, td(R, C, 0), 0))


def k_max_eigenvalues_eigenvectors(k, Matrix):
    lm, v = np.linalg.eigh(Matrix)
    U = v[:, ::-1]
    return U[:, :k]


# expectation value
ad_a = np.kron(a.T @ a, I_c)
ad_ad_a_a = np.kron(a.T @ a.T @ a @ a, I_c)
# parameters values

exp_n = []
exp_e = []
exp_n_2 = []
err_step = []
chi = 60
chi_prime = 4
tau = 0.01
pulse_step = int(0.2 / tau)
Lambda, dw, kappa = 7.5 * np.pi, 2 * np.pi, 0.02 * np.pi
Omega_0 = 1.3 * np.pi
gamma = kappa / 40
drive = 40 * kappa
# Lambda, dw, kappa = 9.5 * np.pi, 2 * np.pi, 0.02 * np.pi
# Omega_0 = 1.5 * np.pi
# gammah = kappa / 400
# gammap = 33 * gammah * 0
# drive = 0
Domega = []
G = []
fn = 7  # 频梳数目
Np = int(N / fn)
mp = int(m / fn)
for n in range(fn):
    for i in range(mp):
        jdw = (n - int(fn / 2)) * dw
        gj = Omega_0 * np.exp(-jdw ** 2 / (2 * Lambda ** 2)) / np.sqrt(Np)
        Domega.append(jdw)
        G.append(gj)

down = np.array(qt.basis(Dms, Dms - 1))
vac = np.array(qt.fock(Dmc, 0))
rho_e = np.kron(down, down)
rho_c = np.kron(vac, vac)
Ide_e = I_e.reshape(-1, 1)
Ide_c = I_c.reshape(-1, 1)
coff = rho_e
Trans_L = [0] * (m - 2)
Trans_R = [0] * (m - 2)
Trace_L = [Ide_e] + [0] * (m - 2)
Trace_R = [0] * (m - 2) + [Ide_e]
# initial process
for n in range(m - 2):
    coff = np.kron(coff, rho_e)
    I_L = np.kron(Trace_L[n], Ide_e)
    dim = coff.shape[0]
    if chi_prime < dim:
        rho = update_coff(coff, rho_e, rho_c).reshape(dim, -1)
        rho_L = rho.dot(rho.T.conj())
        U = k_max_eigenvalues_eigenvectors(chi_prime, rho_L)
        I_L = U.T.conj().dot(I_L)
        coff = U.T.conj().dot(coff)
    else:
        U = np.eye(dim)
    Trans_L[n] = U
    Trace_L[n + 1] = I_L

coff = update_coff(coff, rho_e, rho_c)
Exp_Lindblads_1 = Generate_Exp_Lindblads(0.01, N_p, drive, Domega, G, gamma, kappa)
Exp_Lindblads_2 = Generate_Exp_Lindblads(0.01, N_p, 0, Domega, G, gamma, kappa)

for step in range(0, 1000):
    if step < pulse_step:
        ELs_1 = Exp_Lindblads_1
    else:
        ELs_1 = Exp_Lindblads_2
    err = 0
    coff = coff.reshape(-1, Ds2 * 25)
    # 演化Block_R中第N个原子
    coff = (coff @ ELs_1[-1]).reshape(-1, Ds2, 25)
    # 左扫至第1个原子
    for n in range(m - 2):
        index = m - 2 - n
        dim_R = coff.shape[1]
        # 从Block_L中释放出第N-1个原子
        coff = td(Trans_L[index - 1], coff, [1, 0]).reshape(-1, Ds2, dim_R, 25)
        # 将第N-1个原子与Block_R指标交换
        coff = coff.swapaxes(1, 2)
        # 演化第N-1原子
        coff = (coff.reshape(-1, Ds2 * 25) @ ELs_1[index % m]).reshape(-1, dim_R, Ds2, 25)
        # 再次将第N-1个原子与Block_R指标交换
        coff = coff.swapaxes(1, 2)
        # 将第N-1个原子与Block_R重整化
        coff = coff.reshape(-1, Ds2 * dim_R, 25)
        sr = min(Ds2 * dim_R, chi)
        I_R = np.kron(Ide_e, Trace_R[index])
        if sr != Ds2 * dim_R:
            rho = coff.swapaxes(0, 1).reshape(Ds2 * dim_R, -1)
            rho_R = rho @ rho.T.conj()
            U = k_max_eigenvalues_eigenvectors(sr, rho_R)
            coff = td(U.T.conj(), coff, [1, 1]).swapaxes(0, 1)
            I_R = U.T.conj() @ I_R
            Trans_R[index - 1] = U
        else:
            Trans_R[index - 1] = np.eye(sr)
        Trace_R[index - 1] = I_R
        # 归一化向量化的密度矩阵
        Trace = np.squeeze(td(td(Trace_L[index - 1], Trace_R[index - 1], 0), Ide_c, 0))
        norm = td(coff, Trace.conj(), [[0, 1, 2], [0, 1, 2]])
        coff = coff / norm
        err += (1 - norm)

    # 将第1个原子与Block_R交换并演化第一个原子
    coff = coff.swapaxes(0, 1)
    coff = (coff.reshape(-1, Ds2 * 25) @ ELs_1[0]).reshape(-1, Ds2, 25)
    coff = coff.swapaxes(0, 1)
    # 右扫至第N个原子
    for n in range(m - 2):
        dim_L = coff.shape[0]
        index = n + 1
        # 从Block_R中释放第2个原子
        coff = td(Trans_R[index - 1], coff, [1, 1]).swapaxes(0, 1).reshape(dim_L, Ds2, -1, 25)
        # 将第2个原子与剩下的Block_R指标交换
        coff = coff.swapaxes(1, 2)
        # 演化第2原子
        coff = (coff.reshape(-1, Ds2 * 25) @ ELs_1[index % m]).reshape(dim_L, -1, Ds2, 25)
        # 再次将第2个原子与剩下的Block_R指标交换
        coff = coff.swapaxes(1, 2)
        # 将第2个原子与Block_L重整化
        coff = coff.reshape(dim_L * Ds2, -1, 25)
        sl = min(dim_L * Ds2, chi)
        I_L = np.kron(Trace_L[n], Ide_e)
        if sl != dim_L * Ds2:
            rho = coff.reshape(dim_L * Ds2, -1)
            rho_L = rho @ rho.T.conj()
            U = k_max_eigenvalues_eigenvectors(sl, rho_L)
            coff = td(U.T.conj(), coff, [1, 0])
            I_L = U.T.conj() @ I_L
            Trans_L[index - 1] = U
        else:
            Trans_L[index - 1] = np.eye(sl)
        Trace_L[index] = I_L
        # 归一化向量化的密度矩阵
        Trace = np.squeeze(td(td(Trace_L[index], Trace_R[index], 0), Ide_c, 0))
        norm = td(coff, Trace.conj(), [[0, 1, 2], [0, 1, 2]])
        coff = coff / norm
        err += (1 - norm)

    # 演化Block_R中第N个原子
    coff = coff.reshape(-1, Ds2 * 25)
    coff = (coff @ ELs_1[-1]).reshape(-1, Ds2, 25)
    # 计算期望值
    # Coff = coff.copy()
    # for i in range(2):
    #     dim_R = Coff.shape[1]
    #     Coff = td(Trans_L[-1 - i], Coff, [1, 0]).reshape(-1, Ds2 * dim_R, 25)
    #     Coff = td(Coff, Trans_R[-1 - i].T.conj(), [1, 1]).swapaxes(1, 2)
    # dim_R = Coff.shape[1]
    # Coff = td(Trans_L[-3], Coff, [1, 0]).reshape(-1, Ds2, dim_R, 25).swapaxes(0, 1)
    # Trace_1 = np.squeeze(td(td(td(Trace_L[2], Ide_e, 0), Trace_R[-3], 0), Ide_c, 0)).swapaxes(0, 1)
    # value_e = td(Trace_1.reshape(Ds2, -1).conj(), np.kron(jp @ jm, I_e) @ Coff.reshape(Ds2, -1), [[0, 1], [0, 1]])
    # exp_e.append(value_e)

    Trace = np.squeeze(td(td(Trace_L[-1], Trace_R[-1], 0), Ide_c, 0))
    value = np.einsum('abc, cd, abd', Trace.conj(), ad_a, coff)
    value_1 = np.einsum('abc, cd, abd', Trace.conj(), ad_ad_a_a, coff)
    exp_n.append(value)
    exp_n_2.append(value_1)
    err_step.append(err)
    print('completion rates = {:.1%}'.format(step / 1000))

sio.savemat('Exp_n_N_%s_D_%s_new.mat' % (N,chi), {'exp_n_N_%s_D_%s_new' % (N,chi): exp_n, 'exp_e_N_%s_D_%s_new' % (N,chi): exp_e, 'err_N_%s_D_%s_new' % (N,chi): err_step})
t_f = time.time()
seconds = t_f - t_i
m, s = divmod(seconds, 60)
h, m = divmod(m, 60)
print("%02d:%02d:%02d" % (h, m, s))

t = [0.01 * i for i in range(1, 1001)]
# plt.plot(t, exp_n[:20])
# plt.show()

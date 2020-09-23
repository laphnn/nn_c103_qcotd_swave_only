import os, sys
BMat_path='/Users/walkloud/work/research/c51/x_files/code/pythib'
sys.path.append(BMat_path)
import BMat
import gvar as gv
import numpy as np
import lsqfit
import argparse
import matplotlib.pyplot as plt
import matplotlib
from numpy.linalg import cond
from scipy.optimize import least_squares
from scipy.optimize import curve_fit
from scipy import special as scsp

import read_pickle_results

def run_from_ipython():
    try:
        __IPYTHON__
        return True
    except NameError:
        return False

parser = argparse.ArgumentParser(description='Perform ERE analysis')
parser.add_argument('--nn_type', default='deuteron',
    help='deuteron or dineutron? [%(default)s]')
parser.add_argument('--fit_type',        default='xy_ols',
    help='fit type: y_only, xy_ols, [%(default)s]')
parser.add_argument('--Nbs', type=int,
    help='set number of bs samples to Nbs')
parser.add_argument('--vs_mpi',          default=True,  action='store_const', const=False,
    help='scale qcotd and qsq by mpi? [%(default)s]')
parser.add_argument('--bs_bias_correct', default=True,  action='store_const', const=False,
    help='shift bs mean to b0? [%(default)s]')
parser.add_argument('--all_bs',          default=False, action='store_const', const=True,
    help='plot all bs samples in qcotd curves? [%(default)s]')
parser.add_argument('--nplqcd',          default=False, action='store_const', const=True,
    help='show nplqcd results? [%(default)s]')
args = parser.parse_args()
print(args)

bs_bias_correct = args.bs_bias_correct
all_bs          = args.all_bs
plot_nplqcd     = args.nplqcd
fit_type        = args.fit_type
vs_mpi          = args.vs_mpi
nn_str          = args.nn_type

if nn_str == 'deuteron':
    p_path = 'result/posterior_singlet_t05_td10_N_n2_t_5_20_R_n2_t_5_15_ratio_True.pickle_bs'
elif nn_str == 'dineutron':
    p_path = 'result/posterior_triplet_t05_td10_N_n2_t_5_20_R_n2_t_5_15_ratio_True.pickle_bs'

mpi = gv.gvar('0.310810(95)')
L=48

priors = dict()
priors['mainv'] = gv.gvar(0.1,.3)
priors['r'] = gv.gvar(10,100)
priors['P'] = gv.gvar(1,1000)

momRay = {0 : 'ar', 1 : 'oa', 2 : 'pd', 3 : 'cd', 4 : 'oa'}
irreps_clrs = {
    'T1g':'k', 'A2':'b', 'E':'r', 'B1':'g', 'B2':'magenta',
    'A1g':'k', 'A1':'b',
}
level_mrkr = {0:'s', 1:'o', 2:'d', 3:'p', 4:'h', 5:'8', 6:'v'}

fit_results = gv.load(p_path)
try:
    mN = np.array(fit_results[((('0', 'T1g', 0), 'N', '0'), 'e0')])
except:
    mN = np.array(fit_results[((('0', 'A1g', 0), 'N', '0'), 'e0')])

print('m_N = %s' %gv.gvar(mN[0],mN[1:].std()))

if nn_str == 'deuteron':
    mN = np.array(fit_results[((('0', 'T1g', 0), 'N', '0'), 'e0')])
    def isZero(JtimesTwo, Lp, SptimesTwo, chanp, L, StimesTwo, chan):
        return not (JtimesTwo==2
                and Lp==0 and L==0
                and chanp==0 and chan==0
                and SptimesTwo==2 and StimesTwo==2)
elif nn_str == 'dineutron':
    mN = np.array(fit_results[((('0', 'A1g', 0), 'N', '0'), 'e0')])
    def isZero(JtimesTwo, Lp, SptimesTwo, chanp, L, StimesTwo, chan):
        return not (JtimesTwo==0
                and Lp==0 and L==0
                and chanp==0 and chan==0
                and SptimesTwo==0 and StimesTwo==0)
def calcFunc(self, JtimesTwo, Lp, SptimesTwo, chanp, L, StimesTwo, chan, Ecm_over_mref, pSqFuncList):
    return 0.
Kinv = BMat.KMatrix(calcFunc, isZero)
chanList = [BMat.DecayChannelInfo('n','n',1,1,True,True),]

plt.ion()
plt.figure('qcotd',figsize=(7,4))
ax = plt.axes([0.12,0.12,0.87,0.87])

energies_0, energies_bs, qcotd_0, qcotd_bs, qsq_0, qsq_bs = \
    read_pickle_results.get_data(nn_str, p_path, bs_bias_correct=bs_bias_correct, Mbs=args.Nbs)
Nbs = energies_bs['m_n'].shape[0]

if vs_mpi:
    rescale = (mN[0]/mpi).mean
else:
    rescale = 1.

x_dict = dict()
y_dict = dict()
print('Data used in fit')
for k in qcotd_0:
    Psq,irrep,state = k[0][0]
    Psq = int(Psq)
    k_n = '%d_%s_%d' %(Psq,irrep,state)
    # make qSq data
    x_dict[k_n] = dict()
    x_dict[k_n]['irrep']   = irrep
    x_dict[k_n]['Psq']     = Psq
    x_dict[k_n]['qsq0']    = qsq_0[k]
    x_dict[k_n]['m_N']     = mN[0]
    x_dict[k_n]['L']       = L
    x_dict[k_n]['mom_ray'] = momRay[Psq]

    y_dict[k_n] = qsq_bs[k] - qsq_bs[k].mean() + qsq_0[k]

    clr = irreps_clrs[irrep]
    mkr = level_mrkr[Psq]
    s1,s2 = k[0][2]
    st1 = ((k[0][0], 'N', s1), 'e0')
    st2 = ((k[0][0], 'N', s2), 'e0')
    en1 = np.array(fit_results[st1])
    en2 = np.array(fit_results[st2])
    de_nn0    = energies_0[k]-en1[0]-en2[0]
    de_nnbs   = energies_bs[k]-en1[1:]-en2[1:]
    E_cmSq    = energies_0[k]**2 - Psq*(2*np.pi/L)**2
    E_cmSq_bs = energies_bs[k]**2 - Psq*(2*np.pi/L)**2
    gv_en1    = gv.gvar(en1[0],          en1[1:].std())
    gv_en2    = gv.gvar(en2[0],          en2[1:].std())
    gv_de_nn  = gv.gvar(de_nn0,          de_nnbs.std())
    gv_e_nn   = gv.gvar(energies_0[k],   energies_bs[k].std())
    gv_E_cm   = gv.gvar(np.sqrt(E_cmSq), np.sqrt(E_cmSq_bs).std())
    gv_qsq    = gv.gvar(qsq_0[k],        (qsq_bs[k]/mN[1:]**2).std()) *rescale**2
    gv_qcotd  = gv.gvar(qcotd_0[k],      qcotd_bs[k].std()) *rescale
    print('%d& %3s& %s& %s& %s& %s& %s& %s& %s& %s& %s& %s\\\\' \
        %(Psq, irrep, state, s1, gv_en1, s2, gv_en2, gv_de_nn, gv_e_nn, gv_E_cm, gv_qsq, gv_qcotd))
    i_16 = int(Nbs/100*16)
    i_84 = int(Nbs/100*84)
    q_sort = qcotd_bs[k].argsort()
    ax.plot(qsq_bs[k][q_sort][i_16:i_84]*rescale**2, qcotd_bs[k][q_sort][i_16:i_84]*rescale,
        linestyle='None', color=clr, mfc='None', marker='.', alpha=0.1)
    if state == 0:
        lbl = r'${\rm %s}(P_{\rm  tot}^2 = %d)$' %(irrep,Psq)
    else:
        lbl = ''
    ax.plot(qsq_0[k] * rescale**2, qcotd_0[k] * rescale,
        linestyle='None', color=clr, marker=mkr, label=lbl)

y_gv = gv.dataset.avg_data(y_dict, bstrap=True)

class qsqFit:

    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.irreps = []
        for k in self.x:
            self.irreps.append(k)
            mom_ray = x[k]['mom_ray']
            Psq     = x[k]['Psq']
            irrep   = x[k]['irrep']
            self.x[k]['boxQ'] = BMat.BoxQuantization(mom_ray, Psq, irrep, chanList, [0,], Kinv, True)
            self.x[k]['boxQ'].setRefMassL(x[k]['m_N'] * x[k]['L'])
            self.x[k]['boxQ'].setMassesOverRef(0, 1, 1)

    def ere1(self,x,p0,p1):
        return p0 + 0.5*p1*x
    def ere2(self,x,p0,p1,p2):
        return p0 + 0.5*p1*x + p2*x**2 / 6

    def get_qsq_ere1(self,xl,p0,p1):
        results = []
        for k in xl:
            def residual_sq(x):
                ecm_mn = 2*np.sqrt(1 + x) # x = qSq / mNSq
                qcotd_mn = self.x[k]['boxQ'].getBoxMatrixFromEcm(ecm_mn).real
                res = self.ere1(x,p0,p1) - qcotd_mn
                return res**2

            results.append(least_squares(
                residual_sq, self.x[k]['qsq0'], method='lm', ftol=1.0e-12, gtol=1.0e-12, xtol=1.0e-12
                ).x[0])
        return np.array(results)

    def get_qsq_ere2(self,xl,p0,p1,p2):
        results = []
        for k in xl:
            def residual_sq(x):
                ecm_mn = 2*np.sqrt(1 + x) # x = qSq / mNSq
                qcotd_mn = self.x[k]['boxQ'].getBoxMatrixFromEcm(ecm_mn).real
                res = self.ere2(x,p0,p1,p2) - qcotd_mn
                return res**2

            results.append(least_squares(
                residual_sq, self.x[k]['qsq0'], method='lm', ftol=1.0e-12, gtol=1.0e-12, xtol=1.0e-12
                ).x[0])
        return np.array(results)

    def fit_ere1(self,p0,p1):
        cov_np = np.array(gv.evalcov([self.y[k] for k in self.irreps]))
        yq     = np.array([self.y[k].mean for k in self.irreps])

        p_opt, p_cov = curve_fit(self.get_qsq_ere1, self.irreps, yq, p0=[p0,p1], sigma=cov_np, absolute_sigma=True, method='lm')

        r = self.get_qsq_ere1(self.irreps, p_opt[0],p_opt[1]) - yq
        chisq_min = np.dot(r, np.dot(np.linalg.inv(cov_np), r))
        dof = len(yq) - 2

        results = dict()
        results['chisq_dof'] = chisq_min/dof
        results['dof']       = dof
        results['Q']         = scsp.gammaincc(0.5*dof,0.5*chisq_min)
        results['p_opt']     = gv.gvar(p_opt, p_cov)
        return results

    def fit_ere2(self,p0,p1,p2):
        cov_np = np.array(gv.evalcov([self.y[k] for k in self.irreps]))
        yq     = np.array([self.y[k].mean for k in self.irreps])

        p_opt, p_cov = curve_fit(self.get_qsq_ere2, self.irreps, yq, p0=[p0,p1,p2], sigma=cov_np, absolute_sigma=True, method='lm')

        r = self.get_qsq_ere2(self.irreps, p_opt[0],p_opt[1],p_opt[2]) - yq
        chisq_min = np.dot(r, np.dot(np.linalg.inv(cov_np), r))
        dof = len(yq) - 3

        results = dict()
        results['chisq_dof'] = chisq_min/dof
        results['dof']       = dof
        results['Q']         = scsp.gammaincc(0.5*dof,0.5*chisq_min)
        results['p_opt']     = gv.gvar(p_opt, p_cov)
        return results

qsq_fit = qsqFit(x_dict, y_gv)
p1 = np.array([priors[k].mean for k in ['mainv','r']])
ere1_fit = qsq_fit.fit_ere1(p1[0],p1[1])
print(ere1_fit)
print('\nNLO ERE analysis')
print('  chisq/dof[dof] = %.4f[%d],  Q = %.3f' %(ere1_fit['chisq_dof'],ere1_fit['dof'],ere1_fit['Q']))
mainv = ere1_fit['p_opt'][0]
r     = ere1_fit['p_opt'][1]
print('  -1/am = %s' %(mainv * rescale))
print('     am = %s' %(-1/(mainv * rescale)))
print('   r0 m = %s' %(ere1_fit['p_opt'][1] / rescale))
qp = (1 + np.sqrt(1 + 2 * r * mainv )) / r * rescale
qm = (1 - np.sqrt(1 + 2 * r * mainv )) / r * rescale
print('     qp = %s i' %qp)
print('     qm = %s i' %qm)

p2 = np.array([priors[k].mean for k in ['mainv','r','P']])
ere2_fit = qsq_fit.fit_ere2(p2[0],p2[1],p2[2])
print('\nNNLO ERE analysis')
print('  chisq/dof[dof] = %.4f[%d],  Q = %.3f' %(ere2_fit['chisq_dof'],ere2_fit['dof'],ere2_fit['Q']))
print('  -1/am = %s' %(ere2_fit['p_opt'][0] * rescale))
print('     am = %s' %(-1/(ere2_fit['p_opt'][0] * rescale)))
print('   r0 m = %s' %(ere2_fit['p_opt'][1] / rescale))
print('   r1 m = %s' %(ere2_fit['p_opt'][2] / rescale**3))

if vs_mpi:
    qsq_mN_plot = np.arange(-0.13,0.26,.0005)
else:
    qsq_mN_plot = np.arange(-0.03,0.0605,.0005)

if vs_mpi:
    x_phys = np.arange(-0.13,0.,.00001)
else:
    x_phys = np.arange(-0.04,0.,.00001)
y_phys = -np.sqrt(-x_phys)
ax.plot(x_phys,y_phys,color='cyan')
ax.axhline(color='k')
ax.axvline(color='k')

# q^4 fit
qcotd_m = qsq_fit.ere2(qsq_mN_plot,ere2_fit['p_opt'][0],ere2_fit['p_opt'][1],ere2_fit['p_opt'][2])
y  = np.array([k.mean for k in qcotd_m]) * rescale
dy = np.array([k.sdev for k in qcotd_m]) * rescale
ax.fill_between(qsq_mN_plot*rescale**2,y-dy,y+dy,color='k',alpha=.2)
# q^2 fit
qcotd_m = qsq_fit.ere1(qsq_mN_plot,ere1_fit['p_opt'][0],ere1_fit['p_opt'][1])
y  = np.array([k.mean for k in qcotd_m]) * rescale
dy = np.array([k.sdev for k in qcotd_m]) * rescale
ax.fill_between(qsq_mN_plot*rescale**2,y-dy,y+dy,color='m',alpha=.3)



if vs_mpi:
    ax.set_xlabel(r'$q_{\rm cm}^2 / m_\pi^2$', fontsize=16)
    ax.set_ylabel(r'$q {\rm cot} \delta / m_\pi$', fontsize=16)
else:
    ax.set_xlabel(r'$q_{\rm cm}^2 / m_N^2$', fontsize=16)
    ax.set_ylabel(r'$q {\rm cot} \delta / m_N$', fontsize=16)
if vs_mpi:
    t_cut = (1/2)**2
else:
    t_cut = ((mpi / mN[0] / 2)**2).mean
ax.axvline(t_cut,color='k',linestyle='--')
if vs_mpi:
    ax.axis([-.12, 0.26, -.4,1.2])
else:
    ax.axis([-.026, 0.0525, -.15,0.6])

if not plot_nplqcd:
    ax.legend(loc=2, ncol=5, columnspacing=0, handletextpad=0.1)
else:
    import nplqcd
    clrs = {800:'g'}#,450:'k'} #450 data has issues
    for clr in clrs:
        qsq   = nplqcd.qsq[clr]
        qcotd = nplqcd.qcotd[clr]
        if vs_mpi:
            norm = nplqcd.mpi[clr]
        else:
            norm = nplqcd.mn[clr]
        qsq = qsq / norm**2
        qcotd = qcotd / norm
        y = [k.mean for k in qcotd]
        dy = [k.sdev for k in qcotd]
        x  = [k.mean for k in qsq]
        dx = [k.sdev for k in qsq]
        ax.errorbar(x,y, yerr=dy, xerr=dx,linestyle='None',marker='s', color=clrs[clr])
        ax.annotate(r'$m_\pi\sim%d$: single stout' %clr,
            xy=nplqcd.xy[clr], xycoords='data',
            xytext=nplqcd.xyt[clr], textcoords='data',
            arrowprops=dict(arrowstyle="simple",color=clrs[clr],
                            connectionstyle="arc3,rad=-0.2"),
            fontsize=20,color=clrs[clr]
            )
    ax.annotate(r'$m_\pi\sim714$: CLS',
        xy=(0.08,0.49), xycoords='data',
        xytext=(-0.115, 1.0), textcoords='data',
        arrowprops=dict(arrowstyle="simple",color='magenta',alpha=0.3,
                        connectionstyle="arc3,rad=-0.2"),
        fontsize=20,color='magenta'
        )


plt.figure('qcotd')
if not os.path.exists('plots/qcotd'):
    os.makedirs('plots/qcotd')
if plot_nplqcd:
    plt.savefig(f'plots/qcotd/{nn_str}_qcotd_ere_nplqcd.pdf', transparent=True)
else:
    plt.savefig(f'plots/qcotd/{nn_str}_qcotd_ere_inverse.pdf', transparent=True)


plt.ioff()
if run_from_ipython():
    plt.show(block=False)
else:
    plt.show()

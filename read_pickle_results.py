import numpy as np
import gvar as gv
import os, sys
BMat_path='/Users/walkloud/work/research/c51/x_files/code/pythib'
sys.path.append(BMat_path)
import BMat

L=48
irreps_clrs = {
    'T1g':'k', 'A2':'b', 'E':'r', 'B1':'g', 'B2':'magenta',
    'A1g':'k', 'A1':'b',
}

level_mrkr = {0:'s', 1:'o', 2:'d', 3:'p', 4:'h', 5:'8', 6:'v'}

momRay = {0 : 'ar', 1 : 'oa', 2 : 'pd', 3 : 'cd', 4 : 'oa'}
def calcFunc(self, JtimesTwo, Lp, SptimesTwo, chanp, L, StimesTwo, chan, Ecm_over_mref, pSqFuncList):
    return 0.
chanList = [BMat.DecayChannelInfo('n','n',1,1,True,True),]


def get_data(nn_type,path,bs_bias_correct=True,Mbs=None):
    if nn_type == 'deuteron':
        def isZero(JtimesTwo, Lp, SptimesTwo, chanp, L, StimesTwo, chan):
            return not (JtimesTwo==2
                    and Lp==0 and L==0
                    and chanp==0 and chan==0
                    and SptimesTwo==2 and StimesTwo==2)
    elif nn_type == 'dineutron':
        def isZero(JtimesTwo, Lp, SptimesTwo, chanp, L, StimesTwo, chan):
            return not (JtimesTwo==0
                    and Lp==0 and L==0
                    and chanp==0 and chan==0
                    and SptimesTwo==0 and StimesTwo==0)
    Kinv = BMat.KMatrix(calcFunc, isZero)
    fit_results = gv.load(path)
    try:
        mN = np.array(fit_results[((('0', 'T1g', 0), 'N', '0'), 'e0')])
    except:
        mN = np.array(fit_results[((('0', 'A1g', 0), 'N', '0'), 'e0')])

    n_bs = len(mN[1:])
    if Mbs is not None:
        if Mbs > n_bs:
            sys.exit('you only have '+str(n_bs)+' bs total samples')
        Nbs = Mbs
        mN = mN[:Nbs+1]
    else:
        Nbs = n_bs

    energies_0  = dict()
    energies_bs = dict()
    qcotd_0     = dict()
    qcotd_bs    = dict()
    qsq_0       = dict()
    qsq_bs      = dict()

    energies_0['m_n'] = np.array(mN[0])
    energies_bs['m_n'] = np.array(mN[1:])

    excluded = []

    for k in fit_results:
        if k[1] == 'e0' and k[0][1] == 'R' and len(k[0][0]) == 3 and k[0]:
            Psq, irrep, n = k[0][0]
            Psq = int(Psq)
            if irrep in irreps_clrs and n <= 7:
                if irrep in ['B1','B2'] and Psq != 2:
                    pass
                else:
                    clr = irreps_clrs[irrep]
                    mkr = level_mrkr[Psq]
                    de_nn = np.array(fit_results[k])[:Nbs+1]
                    s1,s2 = k[0][2]
                    st1 = ((k[0][0], 'N', s1), 'e0')
                    st2 = ((k[0][0], 'N', s2), 'e0')
                    en1 = np.array(fit_results[st1])[:Nbs+1]
                    en2 = np.array(fit_results[st2])[:Nbs+1]
                    e_nn = de_nn + en1 + en2
                    E_cmSq = e_nn**2 - Psq*(2*np.pi/L)**2
                    qsq = E_cmSq / 4 - mN**2

                    # make mean (boot0) values
                    mN_0   = mN[0]
                    e_nn_0 = e_nn[0]
                    boxQ = BMat.BoxQuantization(momRay[Psq], Psq, irrep, chanList, [0,], Kinv, True)
                    boxQ.setRefMassL(mN[0]*L)
                    boxQ.setMassesOverRef(0, 1, 1)
                    qcotd = boxQ.getBoxMatrixFromElab(e_nn_0 / mN_0)# in mN units

                    # make bootstrap distribution
                    qsq_qcotd_bs = np.zeros([Nbs,2])
                    boxQ_bs = BMat.BoxQuantization(momRay[Psq], Psq, irrep, chanList, [0,], Kinv, True)
                    boxQ_bs.setMassesOverRef(0,1,1)

                    mN_bs     = mN[1:]
                    e_nn_bs   = e_nn[1:]
                    if bs_bias_correct:
                        dmN     = mN_bs - mN_bs.mean()
                        mN_bs   = mN_0 + dmN
                        denn    = e_nn_bs - e_nn_bs.mean()
                        e_nn_bs = e_nn_0 + denn

                    E_cmSq_bs = e_nn_bs**2 - Psq*(2*np.pi/L)**2
                    qsq_qcotd_bs[:,0] = E_cmSq_bs/4 - mN_bs**2
                    for bs in range(Nbs):
                        boxQ_bs.setRefMassL(mN_bs[bs]*L)
                        qsq_qcotd_bs[bs,1] = boxQ_bs.getBoxMatrixFromElab(e_nn_bs[bs] / mN_bs[bs]).real

                    # save data for ERE analysis
                    gv_en1   = gv.gvar(en1[0],             en1[1:].std())
                    gv_en2   = gv.gvar(en2[0],             en2[1:].std())
                    gv_de_nn = gv.gvar(de_nn[0],           de_nn[1:].std())
                    gv_e_nn  = gv.gvar(e_nn[0],            e_nn[1:].std())
                    gv_E_cm  = gv.gvar(np.sqrt(E_cmSq[0]), np.sqrt(E_cmSq[1:]).std())
                    gv_qsq   = gv.gvar(qsq[0]/mN[0]**2,    (qsq[1:]/mN_bs**2).std())
                    gv_qcotd = gv.gvar(qcotd.real, qsq_qcotd_bs[:,1].std())

                    if qsq[0]/mN_0**2 < 0.05 and np.real(qcotd) < 0.8 and np.real(qcotd) > -0.2:
                        energies_0[k] = e_nn[0]
                        qcotd_0[k]    = qcotd.real
                        qsq_0[k]      = qsq[0]/mN[0]**2
                        #print('%d& %3s& %s& %s& %s& %s& %s& %s& %s& %s& %s& %s\\\\' \
                        #    %(Psq, irrep, n, s1, gv_en1, s2, gv_en2, gv_de_nn, gv_e_nn, gv_E_cm, gv_qsq, gv_qcotd))

                        energies_bs[k] = e_nn_bs
                        qcotd_bs[k]    = qsq_qcotd_bs[:,1]
                        qsq_bs[k]      = qsq_qcotd_bs[:,0] / mN[1:]**2

                    else:
                        excluded.append('%d& %3s& %s& %s& %s& %s& %s& %s& %s& %s& %s& %s\\\\' \
                            %(Psq, irrep, n, s1, gv_en1, s2, gv_en2, gv_de_nn, gv_e_nn, gv_E_cm, gv_qsq, gv_qcotd))

                    '''
                    # make sorted BS values to make plot
                    qsq_qcotd_bs_sorted = qsq_qcotd_bs[qsq_qcotd_bs[:,1].argsort()]

                    qcotd_bs_s = qsq_qcotd_bs_sorted[:,1]
                    qsq_bs_s   = qsq_qcotd_bs_sorted[:,0]
                    if all_bs:
                        i_16=0
                        i_84=-1
                    else:
                        i_16 = int(Nbs/100*16)
                        i_84 = int(Nbs/100*84)

                    ax.plot(qsq_bs_s[i_16:i_84]/mN[0]**2*rescale**2, qcotd_bs_s[i_16:i_84]*rescale,
                        linestyle='None', color=clr, mfc='None', marker='.', alpha=0.1)

                    if n == 0:
                        lbl = r'${\rm %s}(P_{\rm  tot}^2 = %d)$' %(irrep,Psq)
                    else:
                        lbl = ''
                    ax.plot(qsq[0]/mN[0]**2 * rescale**2, qcotd.real * rescale,
                        linestyle='None', color=clr, marker=mkr, label=lbl)
                    '''
    print('\nExcluded')
    if len(excluded) == 0:
        print('None')
    else:
        for e in excluded:
            print(e)

    return energies_0, energies_bs, qcotd_0, qcotd_bs, qsq_0, qsq_bs
